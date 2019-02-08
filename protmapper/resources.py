from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import boto3
import logging
import requests
import botocore
from urllib.request import urlretrieve


logger = logging.getLogger(__name__)


# If the protmapper resource directory does not exist, try to create it
home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.protmapper')


if not os.path.isdir(resource_dir):
    try:
        os.makedirs(resource_dir)
    except Exception:
        logger.warning(resource_dir + ' already exists')


def download_phosphositeplus(out_file, cached=True):
    psp_url = ('http://sorger.med.harvard.edu/data/bachman/'
                       'Phosphorylation_site_dataset.tsv')
    logger.info("Downloading PhosphoSitePlus data from %s\n" % psp_url)
    logger.info("Note that PhosphoSitePlus data is not available for "
                "commercial use; please see full terms and conditions at: "
                "https://www.psp.org/staticDownloads")
    resp = requests.get(psp_url)
    # Check the status code
    if resp.status_code == 200:
        # Read and write as bytes (response.content)
        logger.info("Saving PhosphoSitePlus data to %s" % out_file)
        with open(out_file, 'wb') as f:
            f.write(resp.content)
    else:
        logger.error("Error %s occurred downloading PhosphoSitePlus data" %
                     resp.status_code)


def download_uniprot_entries(out_file, cached=True):
    logger.info('Downloading UniProt entries')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'format=tab&columns=id,genes(PREFERRED),' + \
        'entry%20name,database(RGD),database(MGI),length'
    logger.info('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        logger.info('Failed to download "%s"' % url)
    reviewed_entries = res.content

    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:no&fil=organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&' + \
        'format=tab&columns=id,genes(PREFERRED),entry%20name,' + \
        'database(RGD),database(MGI),length'
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
    if not cached:
        logger.info('Downloading UniProt secondary accession mappings')
        url = 'ftp://ftp.uniprot.org/pub/databases/uniprot/knowledgebase/' + \
            'docs/sec_ac.txt'
        urlretrieve(url, out_file)
    else:
        s3 = boto3.resource('s3')
        s3.Bucket('bigmech').download_file('travis/uniprot_sec_ac.txt',
                                           out_file)


def download_hgnc_entries(out_file, cached=True):
    logger.info('Downloading HGNC entries')
    url = 'http://tinyurl.com/y83dx5s6'
    res = requests.get(url)
    if res.status_code != 200:
        logger.error('Failed to download "%s"' % url)
        return
    logger.info('Saving into %s' % out_file)
    with open(out_file, 'wb') as fh:
        fh.write(res.content)


RESOURCE_MAP = {
    'hgnc': ('hgnc_entries.tsv', download_hgnc_entries),
    'upsec': ('uniprot_sec_ac.txt', download_uniprot_sec_ac),
    'up': ('uniprot_entries.tsv', download_uniprot_entries),
    'psp': ('Phosphorylation_site_dataset.tsv', download_phosphositeplus),
    }


class ResourceManager(object):
    def __init__(self, resource_map):
        self.resource_map = resource_map

    def get_resource_file(self, resource_id):
        return os.path.join(resource_dir, self.resource_map[resource_id][0])

    def get_download_fun(self, resource_id):
        return self.resource_map[resource_id][1]

    def has_resource_file(self, resource_id):
        fname = self.get_resource_file(resource_id)
        return os.path.exists(fname)

    def download_resource_file(self, resource_id, cached=True):
        download_fun = self.get_download_fun(resource_id)
        fname = self.get_resource_file(resource_id)
        logger.info('Downloading \'%s\' resource file into %s%s.' %
                    (resource_id, fname, ' from cache' if cached else ''))
        download_fun(fname, cached=cached)

    def get_create_resource_file(self, resource_id, cached=True):
        if not self.has_resource_file(resource_id):
            self.download_resource_file(resource_id, cached)
        return self.get_resource_file(resource_id)

    def get_resource_ids(self):
        return list(self.resource_map.keys())


resource_manager = ResourceManager(RESOURCE_MAP)


if __name__ == '__main__':
    resource_ids = resource_manager.get_resource_ids()
    for resource_id in resource_ids:
        resource_manager.get_create_resource_file(resource_id)
