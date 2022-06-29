import gzip

import os
import re
import csv
import zlib
import json
import boto3
import pystow
import logging
import argparse
import requests
import botocore
from ftplib import FTP
from io import BytesIO, StringIO
from collections import namedtuple
from urllib.request import urlretrieve
from xml.etree import ElementTree as ET
from . import __version__


logger = logging.getLogger('protmapper.resources')
logger.setLevel(logging.INFO)


# If the protmapper resource directory does not exist, try to create it using PyStow
# Can be specified with PROTMAPPER_HOME environment variable, otherwise defaults
# to $HOME/.data/protmapper/<__version__>. The location of $HOME can be overridden with
# the PYSTOW_HOME environment variable
resource_dir_path = pystow.join('protmapper', __version__)
resource_dir = resource_dir_path.as_posix()


def _download_from_s3(key, out_file):
    s3 = boto3.client('s3',
                      config=botocore.client.Config(
                          signature_version=botocore.UNSIGNED))
    tc = boto3.s3.transfer.TransferConfig(use_threads=False)
    # Path to the versioned resource file
    full_key = 'protmapper/%s/%s' % (__version__, key)
    s3.download_file('bigmech', full_key, out_file, Config=tc)


def _download_ftp(ftp_host, ftp_path, out_file=None, ftp_blocksize=33554432,
                  decompress=True):
    ftp = FTP(ftp_host)
    ftp.login()
    fbytes = BytesIO()
    ftp.retrbinary('RETR %s' % ftp_path,
                   callback=lambda s: fbytes.write(s),
                   blocksize=ftp_blocksize)
    ret = fbytes.getvalue()
    if decompress:
        ret = zlib.decompress(ret, 16+zlib.MAX_WBITS)
    if out_file is not None:
        with open(out_file, 'wb') as f:
            f.write(ret)
    return ret


def download_phosphositeplus(out_file, cached=True):
    if not cached:
        logger.info('Cannot download PSP data without using the cache.')
        return
    logger.info("Note that PhosphoSitePlus data is not available for "
                "commercial use; please see full terms and conditions at: "
                "https://www.psp.org/staticDownloads")
    _download_from_s3('Phosphorylation_site_dataset.tsv.gz', out_file)


def download_uniprot_entries(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_entries.tsv.gz', out_file)
        return
    base_columns = ['id', 'genes(PREFERRED)', 'entry%20name',
                    'database(RGD)', 'database(MGI)', 'length', 'reviewed',
                    'organism-id', 'database(GeneID)']
    processed_columns = ['genes', 'protein%20names']
    feature_types = ['SIGNAL', 'CHAIN', 'PROPEPTIDE', 'PEPTIDE', 'TRANSIT']
    columns = base_columns + processed_columns + \
        ['feature(%s)' % feat for feat in feature_types]
    columns_str = ','.join(columns)
    logger.info('Downloading UniProt entries')
    url = 'https://legacy.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'format=tab&columns=' + columns_str
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.info('Failed to download "%s"' % url)
    reviewed_entries = res.content

    url = 'https://legacy.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:no&fil=organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&' + \
        'format=tab&columns=' + columns_str
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.info('Failed to download "%s"' % url)
    unreviewed_human_entries = res.content

    if not((reviewed_entries is not None) and
            (unreviewed_human_entries is not None)):
        return
    unreviewed_human_entries = unreviewed_human_entries.decode('utf-8')
    reviewed_entries = reviewed_entries.decode('utf-8')
    lines = reviewed_entries.strip('\n').split('\n')
    lines += unreviewed_human_entries.strip('\n').split('\n')[1:]

    logger.info('Processing UniProt entries list.')
    new_lines = ['\t'.join(base_columns + processed_columns + ['features'])]
    for line_idx, line in enumerate(lines):
        if line_idx == 0:
            continue
        new_line = process_uniprot_line(line, base_columns, processed_columns,
                                        feature_types)
        new_lines.append(new_line)

    # Join all lines into a single string
    full_table = '\n'.join(new_lines)
    logger.info('Saving into %s.' % out_file)
    with gzip.open(out_file, 'wt', encoding='utf-8') as fh:
        fh.write(full_table)


def process_uniprot_line(line, base_columns, processed_columns,
                         feature_types):
    terms = line.split('\t')

    # At this point, we need to clean up the gene names.
    # If there are multiple gene names, take the first one
    gene_names_preferred = terms[1].split(';')
    gene_name = gene_names_preferred[0]
    if not gene_name:
        gene_name = terms[len(base_columns)].split(' ')[0]

    protein_names = parse_uniprot_synonyms(terms[len(base_columns)+1])
    protein_name = protein_names[0] if protein_names else None

    if gene_name:
        terms[1] = gene_name
    elif protein_name:
        terms[1] = protein_name
    else:
        terms[1] = None

    # We only add Entrez IDs for reviewed entries to avoid the problem
    # caused by one-to-many mappings with lots of unreviewed proteins
    terms[8] = '' if terms[6] != 'reviewed' else terms[8]

    # Next we process the various features into a form that can be
    # loaded easily in the client
    features = []
    for idx, feature_type in enumerate(feature_types):
        col_idx = len(base_columns) + len(processed_columns) + idx
        features += _process_feature(feature_type, terms[col_idx],
                                     protein_name)
    features_json = [feature_to_json(feature) for feature in features]
    features_json_str = json.dumps(features_json)
    new_line = terms[:len(base_columns)] + [features_json_str]
    return '\t'.join(new_line)


def parse_uniprot_synonyms(synonyms_str):
    synonyms_str = re.sub(r'\[Includes: ([^]])+\]',
                          '', synonyms_str).strip()
    synonyms_str = re.sub(r'\[Cleaved into: ([^]])+\]( \(Fragments\))?',
                          '', synonyms_str).strip()

    def find_block_from_right(s):
        parentheses_depth = 0
        assert s.endswith(')')
        s = s[:-1]
        block = ''
        for c in s[::-1]:
            if c == ')':
                parentheses_depth += 1
            elif c == '(':
                if parentheses_depth > 0:
                    parentheses_depth -= 1
                else:
                    return block
            block = c + block
        return block

    syns = []
    while True:
        if not synonyms_str:
            return syns
        if not synonyms_str.endswith(')'):
            return [synonyms_str] + syns

        syn = find_block_from_right(synonyms_str)
        synonyms_str = synonyms_str[:-len(syn)-3]
        # EC codes are not valid synonyms
        if not re.match(r'EC [\d\.-]+', syn):
            syns = [syn] + syns


Feature = namedtuple('Feature', ['type', 'begin', 'end', 'name', 'id',
                                 'is_main'])


def feature_to_json(feature):
    jj = {
        'type': feature.type,
        'begin': feature.begin,
        'end': feature.end,
        'name': feature.name,
        'id': feature.id
    }
    if feature.is_main:
        jj['is_main'] = True
    return jj


def feature_from_json(feature_json):
    if 'is_main' not in feature_json:
        feature_json['is_main'] = False
    return Feature(**feature_json)


def _process_feature(feature_type, feature_str, protein_name):
    """Process a feature string from the UniProt TSV.
    Documentation at: https://www.uniprot.org/help/sequence_annotation
    """
    # This function merges parts that were split inadvertently on semicolons
    def _fix_parts(parts):
        for idx, part in enumerate(parts):
            if part.startswith('/') and not part.endswith('"'):
                parts[idx] += parts[idx+1]
                parts = [p for idx, p in enumerate(parts) if idx != idx+1]
        return parts

    # Split parts and strip off extra spaces
    parts = [p.strip(' ;') for p in feature_str.split('; ')]
    parts = _fix_parts(parts)
    # Find each starting part e.g., CHAIN
    chunk_ids = [idx for idx, part in enumerate(parts)
                 if part.startswith(feature_type)] + [len(parts)]
    # Group parts into chunks, one for each overall entry
    chunks = []
    for idx, chunk_id in enumerate(chunk_ids[:-1]):
        chunks.append(parts[chunk_ids[idx]:chunk_ids[idx+1]])
    feats = []
    # For each distinct entry, we collect all the relevant parts and parse
    # out information
    for chunk in chunks:
        begin = end = name = pid = None
        for part in chunk:
            # If this is the starting piece, we have to parse out the begin
            # and end coordinates. Caveats include: sometimes only one
            # number is given; sometimes a ? is there instead of a number;
            # sometimes a question mark precedes a number; sometimes
            # there is a < before the beginning number; sometimes there
            # is a > before the end number. Sometimes there is an isoform
            # before the beginning number.
            if part.startswith(feature_type):
                match = re.match(r'%s '                           # type marker
                                 r'(?:[^:]+:)?(?:\?|<?)(\d+|\?)'  # beginning
                                 r'..'                            # connector
                                 r'(?:\?|>?)(\d+|\?)' %           # end
                                 feature_type, part)
                if match:
                    beg, end = match.groups()
                else:
                    # This is the standard begin marker
                    match = re.match(r'%s (\d+)' % feature_type, part)
                    beg = match.groups()[0]
                    end = beg
                begin = int(beg) if beg != '?' else None
                end = int(end) if end != '?' else None
            elif part.startswith('/note'):
                match = re.match(r'/note="(.+)"', part)
                name = match.groups()[0]
            elif part.startswith('/id'):
                match = re.match(r'/id="(.+)"', part)
                pid = match.groups()[0]
        is_main = (name == protein_name)
        feature = Feature(feature_type, begin, end, name, pid, is_main)
        feats.append(feature)
    return feats


def download_uniprot_sec_ac(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sec_ac.txt.gz', out_file)
        return

    logger.info('Downloading UniProt secondary accession mappings')
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
        'complete/docs/sec_ac.txt'
    content = _download_ftp('ftp.uniprot.org',
                            'pub/databases/uniprot/knowledgebase/complete/'
                            'docs/sec_ac.txt',
                            decompress=False, out_file=None)
    with gzip.open(out_file, 'wb') as fh:
        fh.write(content)


def download_hgnc_entries(out_file, cached=True):
    if cached:
        _download_from_s3('hgnc_entries.tsv.gz', out_file)
        return

    logger.info('Downloading HGNC entries')
    url = 'http://tinyurl.com/y83dx5s6'
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download "%s"' % url)
        return
    logger.info('Saving into %s' % out_file)
    with gzip.open(out_file, 'wt', encoding='utf-8') as fh:
        fh.write(res.text)


def download_swissprot(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sprot.fasta.gz', out_file)
        return
    logger.info('Downloading reviewed protein sequences from SwissProt')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                'complete/uniprot_sprot.fasta.gz')
    _download_ftp('ftp.uniprot.org', ftp_path, out_file, decompress=False)


def download_isoforms(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sprot_varsplic.fasta.gz', out_file)
        return
    logger.info('Downloading isoform sequences from Uniprot')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                'complete/uniprot_sprot_varsplic.fasta.gz')
    _download_ftp('ftp.uniprot.org', ftp_path, out_file, decompress=False)


def download_refseq_seq(out_file, cached=True):
    if cached:
        _download_from_s3('refseq_sequence.fasta.gz', out_file)
        return
    ftp_path = ('/refseq/H_sapiens/annotation/GRCh38_latest/'
                'refseq_identifiers/GRCh38_latest_protein.faa.gz')
    _download_ftp('ftp.ncbi.nlm.nih.gov', ftp_path, out_file,
                  decompress=False)


def download_refseq_uniprot(out_file, cached=True):
    if cached:
        _download_from_s3('refseq_uniprot.csv.gz', out_file)
        return
    logger.info('Downloading RefSeq->Uniprot mappings from Uniprot')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                'idmapping/by_organism/HUMAN_9606_idmapping.dat.gz')
    mappings_bytes = _download_ftp('ftp.uniprot.org', ftp_path,
                                   out_file=None, decompress=True)
    logger.info('Processing RefSeq->Uniprot mappings file')
    mappings_io = StringIO(mappings_bytes.decode('utf8'))
    csvreader = csv.reader(mappings_io, delimiter='\t')
    filt_rows = []
    for up_id, other_type, other_id in csvreader:
        if other_type == 'RefSeq':
            filt_rows.append([other_id, up_id])
    # Write the file with just the RefSeq->UP mappings
    with gzip.open(out_file, 'wt', encoding='utf-8') as fh:
        writer = csv.writer(fh)
        writer.writerows(filt_rows)


RESOURCE_MAP = {
    'hgnc': ('hgnc_entries.tsv.gz', download_hgnc_entries),
    'upsec': ('uniprot_sec_ac.txt.gz', download_uniprot_sec_ac),
    'up': ('uniprot_entries.tsv.gz', download_uniprot_entries),
    'psp': ('Phosphorylation_site_dataset.tsv.gz', download_phosphositeplus),
    'swissprot': ('uniprot_sprot.fasta.gz', download_swissprot),
    'isoforms': ('uniprot_sprot_varsplic.fasta.gz', download_isoforms),
    'refseq_uniprot': ('refseq_uniprot.csv.gz', download_refseq_uniprot),
    'refseq_seq': ('refseq_sequence.fasta.gz', download_refseq_seq),
    }


class ResourceManager(object):
    """Class to manage a set of resource files.

    Parameters
    ----------
    resource_map : dict
        A dict that maps resource file IDs to a tuple of resource file names
        and download functions.
    """
    def __init__(self, resource_map):
        self.resource_map = resource_map

    def get_resource_file(self, resource_id):
        """Return the path to the resource file with the given ID.

        Parameters
        ----------
        resource_id : str
            The ID of the resource.

        Returns
        -------
        str
            The path to the resource file.
        """
        return os.path.join(resource_dir_path, self.resource_map[resource_id][0])

    def get_download_fun(self, resource_id):
        """Return the download function for the given resource.

        Parameters
        ----------
        resource_id : str
            The ID of the resource.

        Returns
        -------
        function
            The download function for the given resource.
        """
        return self.resource_map[resource_id][1]

    def has_resource_file(self, resource_id):
        """Return True if the resource file exists for the given ID.

        Parameters
        ----------
        resource_id : str
            The ID of the resource.

        Returns
        -------
        bool
            True if the resource file exists, false otherwise.
        """
        fname = self.get_resource_file(resource_id)
        return os.path.exists(fname)

    def download_resource_file(self, resource_id, cached=True):
        """Download the resource file corresponding to the given ID.

        Parameters
        ----------
        resource_id : str
            The ID of the resource.
        cached : Optional[bool]
            If True, the download is a pre-processed file from S3, otherwise
            the download is obtained and processed from the primary source.
            Default: True
        """
        download_fun = self.get_download_fun(resource_id)
        fname = self.get_resource_file(resource_id)
        logger.info('Downloading \'%s\' resource file into %s%s.' %
                    (resource_id, fname, ' from cache' if cached else ''))
        download_fun(fname, cached=cached)

    def get_create_resource_file(self, resource_id, cached=True):
        """Return the path to the resource file, download if it doesn't exist.

        Parameters
        ----------
        resource_id : str
            The ID of the resource.
        cached : Optional[bool]
            If True, the download is a pre-processed file from S3, otherwise
            the download is obtained and processed from the primary source.
            Default: True

        Returns
        -------
        str
            The path to the resource file.
        """
        if not self.has_resource_file(resource_id):
            logger.info(('Could not access \'%s\' resource'
                         ' file, will download.') % resource_id)
            self.download_resource_file(resource_id, cached)
        return self.get_resource_file(resource_id)

    def get_resource_ids(self):
        """Return a list of all the resource IDs managed by this manager."""
        return list(self.resource_map.keys())


resource_manager = ResourceManager(RESOURCE_MAP)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # By default we use the cache
    parser.add_argument('--uncached', action='store_true')
    # By default we use get_create which doesn't do anything if the resource
    # already exists. With the download flag, we force re-download.
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    logger.info(args)
    resource_ids = resource_manager.get_resource_ids()
    for resource_id in resource_ids:
        if not args.download:
            resource_manager.get_create_resource_file(resource_id,
                                                      cached=(not
                                                              args.uncached))
        else:
            resource_manager.download_resource_file(resource_id,
                                                    cached=(not args.uncached))

