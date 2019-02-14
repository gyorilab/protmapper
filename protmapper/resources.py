import os
import boto3
import logging
import requests
import botocore
from urllib.request import urlretrieve


logger = logging.getLogger('protmapper.resources')


# If the protmapper resource directory does not exist, try to create it
home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.protmapper')


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
    s3.download_file('bigmech', 'travis/%s' % key, out_file, Config=tc)


def download_phosphositeplus(out_file, cached=True):
    logger.info("Note that PhosphoSitePlus data is not available for "
                "commercial use; please see full terms and conditions at: "
                "https://www.psp.org/staticDownloads")
    _download_from_s3('Phosphorylation_site_dataset.tsv', out_file)


def download_uniprot_entries(out_file, cached=True):
    if cached:
        _download_from_s3('uniprot_entries.tsv', out_file)
        return

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


RESOURCE_MAP = {
    'hgnc': ('hgnc_entries.tsv', download_hgnc_entries),
    'upsec': ('uniprot_sec_ac.txt', download_uniprot_sec_ac),
    'up': ('uniprot_entries.tsv', download_uniprot_entries),
    'psp': ('Phosphorylation_site_dataset.tsv', download_phosphositeplus),
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
        """Return True if the resource file existis for the given ID.

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
    resource_ids = resource_manager.get_resource_ids()
    for resource_id in resource_ids:
        resource_manager.get_create_resource_file(resource_id)
