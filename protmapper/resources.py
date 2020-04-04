import os
import csv
import zlib
import boto3
import logging
import argparse
import requests
import botocore
from ftplib import FTP
from io import BytesIO, StringIO
from urllib.request import urlretrieve
from xml.etree import ElementTree as ET
from . import __version__


logger = logging.getLogger('protmapper.resources')


# If the protmapper resource directory does not exist, try to create it
home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.protmapper', __version__)


if not os.path.isdir(resource_dir):
    try:
        os.makedirs(resource_dir)
    except Exception:
        logger.warning(resource_dir + ' already exists')


def _download_from_s3(key, out_file):
    s3 = boto3.client('s3',
                      config=botocore.client.Config(
                          signature_version=botocore.UNSIGNED))
    tc = boto3.s3.transfer.TransferConfig(use_threads=False)
    # Path to the versioned resource file
    full_key = 'protmapper/%s/%s' % (__version__, key)
    s3.download_file('bigmech', full_key, out_file, Config=tc)


def _download_ftp_gz(ftp_host, ftp_path, out_file=None, ftp_blocksize=33554432):
    ftp = FTP(ftp_host)
    ftp.login()
    gzf_bytes = BytesIO()
    ftp.retrbinary('RETR %s' % ftp_path,
                   callback=lambda s: gzf_bytes.write(s),
                   blocksize=ftp_blocksize)
    ret = gzf_bytes.getvalue()
    ret = zlib.decompress(ret, 16+zlib.MAX_WBITS)
    if out_file is not None:
        with open(out_file, 'wb') as f:
            f.write(ret)
    return ret


def download_phosphositeplus(out_file, cached=True):
    logger.info("Note that PhosphoSitePlus data is not available for "
                "commercial use; please see full terms and conditions at: "
                "https://www.psp.org/staticDownloads")
    _download_from_s3('Phosphorylation_site_dataset.tsv', out_file)


def download_uniprot_entries(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_entries.tsv', out_file)
        return

    columns = ['id', 'genes(PREFERRED)', 'entry%20name', 'database(RGD)',
               'database(MGI)', 'length', 'reviewed', 'feature(SIGNAL)',
               'feature(CHAIN)', 'feature(PROPEPTIDE)']
    columns_str = ','.join(columns)

    logger.info('Downloading UniProt entries')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'format=tab&columns=' + columns_str
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.info('Failed to download "%s"' % url)
    reviewed_entries = res.content

    url = 'http://www.uniprot.org/uniprot/?' + \
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
    # At this point, we need to clean up the gene names.
    logger.info('Processing UniProt entries list.')
    for i, line in enumerate(lines):
        if i == 0:
            continue
        terms = line.split('\t')
        # If there are multiple gene names, take the first one
        gene_names = terms[1].split(';')
        terms[1] = gene_names[0]
        # Join the line again after the change
        lines[i] = '\t'.join(terms)
    # Join all lines into a single string
    full_table = '\n'.join(lines)
    logging.info('Saving into %s.' % out_file)
    with open(out_file, 'wb') as fh:
        fh.write(full_table.encode('utf-8'))


def download_uniprot_sec_ac(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sec_ac.txt', out_file)
        return

    logger.info('Downloading UniProt secondary accession mappings')
    url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
        'docs/sec_ac.txt'
    urlretrieve(url, out_file)


def download_hgnc_entries(out_file, cached=True):
    if cached:
        _download_from_s3('hgnc_entries.tsv', out_file)
        return

    logger.info('Downloading HGNC entries')
    url = 'http://tinyurl.com/y83dx5s6'
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download "%s"' % url)
        return
    logger.info('Saving into %s' % out_file)
    with open(out_file, 'wb') as fh:
        fh.write(res.content)


def download_swissprot(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sprot.fasta', out_file)
        return
    logger.info('Downloading reviewed protein sequences from SwissProt')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                 'complete/uniprot_sprot.fasta.gz')
    _download_ftp_gz('ftp.uniprot.org', ftp_path, out_file)


def download_isoforms(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sprot_varsplic.fasta', out_file)
        return
    logger.info('Downloading isoform sequences from Uniprot')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                'complete/uniprot_sprot_varsplic.fasta.gz')
    _download_ftp_gz('ftp.uniprot.org', ftp_path, out_file)


def download_refseq_seq(out_file, cached=True):
    if cached:
        _download_from_s3('refseq_sequence.fasta', out_file)
        return
    ftp_path = ('/refseq/H_sapiens/annotation/GRCh38_latest/'
                'refseq_identifiers/GRCh38_latest_protein.faa.gz')
    _download_ftp_gz('ftp.ncbi.nlm.nih.gov', ftp_path, out_file)


def download_refseq_uniprot(out_file, cached=True):
    if cached:
        _download_from_s3('refseq_uniprot.csv', out_file)
        return
    logger.info('Downloading RefSeq->Uniprot mappings from Uniprot')
    ftp_path = ('/pub/databases/uniprot/current_release/knowledgebase/'
                'idmapping/by_organism/HUMAN_9606_idmapping.dat.gz')
    mappings_bytes = _download_ftp_gz('ftp.uniprot.org', ftp_path,
                                      out_file=None)
    logger.info('Processing RefSeq->Uniprot mappings file')
    mappings_io = StringIO(mappings_bytes.decode('utf8'))
    csvreader = csv.reader(mappings_io, delimiter='\t')
    filt_rows = []
    for up_id, other_type, other_id in csvreader:
        if other_type == 'RefSeq':
            filt_rows.append([other_id, up_id])
    # Write the file with just the RefSeq->UP mappings
    with open(out_file, 'wt') as f:
        csvwriter = csv.writer(f)
        csvwriter.writerows(filt_rows)


def download_sars_cov2(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_sars_cov2_entries.tsv', out_file)
    else:
        logger.info('Downloading Sars-Cov-2 mappings from Uniprot')
        url = ('ftp://ftp.uniprot.org/pub/databases/uniprot/pre_release/'
               'covid-19.xml')
        urlretrieve(url, out_file + '.xml')

    et = ET.parse(out_file + '.xml')
    up_ns = {'up': 'http://uniprot.org/uniprot'}

    def _get_chains(entry):
        features = entry.findall('up:feature', namespaces=up_ns)
        chains = []
        for feature in features:
            if feature.attrib.get('type') == 'chain':
                pid = feature.attrib['id']
                desc = feature.attrib['description']
                begin = feature.find('up:location/up:begin',
                                     namespaces=up_ns).attrib['position']
                end = feature.find('up:location/up:end',
                                   namespaces=up_ns).attrib['position']
                chains.append((pid, desc, begin, end))
        return chains

    def make_chain_str(chains):
        chain_strs = ['CHAIN %s..%s;  /note="%s";  /id="%s";' %
                      (begin, end, desc, pid) for pid, desc, begin, end
                      in chains]
        chain_str = '  '.join(chain_strs)
        return chain_str

    rows = [('Entry', 'Gene names  (primary )', 'Entry name',
             'Cross-reference (RGD)', 'Cross-reference (MGI)', 'Length',
             'Status', 'Signal peptide', 'Chain', 'Propeptide')]
    for entry in et.findall('up:entry', namespaces=up_ns):
        up_id = entry.find('up:accession', namespaces=up_ns).text
        mnemonic = entry.find('up:name', namespaces=up_ns).text
        # Skip redundant human proteins here
        if mnemonic.endswith('HUMAN'):
            continue
        gene_name_tag = entry.find('up:gene/up:name', namespaces=up_ns)
        gene_name = gene_name_tag.text if gene_name_tag is not None else None
        full_name_tag = entry.find('up:protein/up:recommendedName/up:fullName',
                                   namespaces=up_ns)
        full_name = full_name_tag.text if full_name_tag is not None else None
        recommended_name_tag = \
            entry.find('up:protein/up:recommendedName/up:fullName',
                       namespaces=up_ns)
        recommended_name = recommended_name_tag.text if recommended_name_tag \
            is not None else None
        short_names = [e.text for e in
                       entry.findall('up:protein/up:recommendedName/'
                                     'up:shortName', namespaces=up_ns)]

        # Choose a single canonical name
        if recommended_name:
            canonical_name = recommended_name
        elif gene_name:
            canonical_name = gene_name
        elif full_name:
            canonical_name = full_name
        elif short_names:
            canonical_name = short_names[0]
        if not canonical_name:
            assert False

        chains = _get_chains(entry)
        chain_str = make_chain_str(chains)
        seq_tag = entry.find('up:sequence', namespaces=up_ns)
        length = seq_tag.attrib['length']

        row = up_id, canonical_name, mnemonic, '', '', length, \
            'reviewed', '', chain_str, ''
        rows.append(row)
    with open(out_file, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', quotechar=None)
        for row in rows:
            writer.writerow(row)


RESOURCE_MAP = {
    'hgnc': ('hgnc_entries.tsv', download_hgnc_entries),
    'upsec': ('uniprot_sec_ac.txt', download_uniprot_sec_ac),
    'up': ('uniprot_entries.tsv', download_uniprot_entries),
    'up_sars_cov2': ('uniprot_sars_cov2_entries.tsv', download_sars_cov2),
    'psp': ('Phosphorylation_site_dataset.tsv', download_phosphositeplus),
    'swissprot': ('uniprot_sprot.fasta', download_swissprot),
    'isoforms': ('uniprot_sprot_varsplic.fasta', download_isoforms),
    'refseq_uniprot': ('refseq_uniprot.csv', download_refseq_uniprot),
    'refseq_seq': ('refseq_sequence.fasta', download_refseq_seq),
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
        return os.path.join(resource_dir, self.resource_map[resource_id][0])

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
    resource_ids = resource_manager.get_resource_ids()
    for resource_id in resource_ids:
        if not args.download:
            resource_manager.get_create_resource_file(resource_id,
                                                      cached=(not
                                                              args.uncached))
        else:
            resource_manager.download_resource_file(resource_id,
                                                    cached=(not args.uncached))

