from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging
import requests


logger = logging.getLogger(__name__)


# If the protmapper resource directory does not exist, try to create it
home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.protmapper')


if not os.path.isdir(resource_dir):
    try:
        os.makedirs(resource_dir)
    except Exception:
        logger.warning(resource + ' already exists')


psp_filename = os.path.join(resource_dir, 'Phosphorylation_site_dataset.tsv')


def download_phosphositeplus():
    psp_url = ('http://sorger.med.harvard.edu/data/bachman/'
                       'Phosphorylation_site_dataset.tsv')
    print("Downloading PhosphoSitePlus data from %s\n" % psp_url)
    print("Note that PhosphoSitePlus data is not available for commercial use; "
          "please see full terms and conditions at: "
          "https://www.psp.org/staticDownloads")
    resp = requests.get(psp_url)
    # Check the status code
    if resp.status_code == 200:
        # Read and write as bytes (response.content)
        logger.info("Saving PhosphoSitePlus data to %s" % psp_filename)
        with open(psp_filename, 'wb') as f:
            f.write(resp.content)
    else:
        logger.error("Error %s occurred downloading PhosphoSitePlus data" %
                     resp.status_code)


def download_uniprot_mappings():
    print('Downloading UniProt entries')
    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:yes&' + \
        'format=tab&columns=id,genes(PREFERRED),' + \
        'entry%20name,database(RGD),database(MGI)'
    print('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        print('Failed to download "%s"' % url)
    reviewed_entries = res.content

    url = 'http://www.uniprot.org/uniprot/?' + \
        'sort=id&desc=no&compress=no&query=reviewed:no&fil=organism:' + \
        '%22Homo%20sapiens%20(Human)%20[9606]%22&' + \
        'format=tab&columns=id,genes(PREFERRED),entry%20name,' + \
        'database(RGD),database(MGI)'
    print('Downloading %s' % url)
    res = requests.get(url)
    if res.status_code != 200:
        print('Failed to download "%s"' % url)
    unreviewed_human_entries = res.content

    if not((reviewed_entries is not None) and
            (unreviewed_human_entries is not None)):
            return
    unreviewed_human_entries = unreviewed_human_entries.decode('utf-8')
    reviewed_entries = reviewed_entries.decode('utf-8')
    lines = reviewed_entries.strip('\n').split('\n')
    lines += unreviewed_human_entries.strip('\n').split('\n')[1:]
    # At this point, we need to clean up the gene names.
    print('Processing UniProt entries list.')
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
    #fname = os.path.join(path, 'uniprot_entries.tsv')
    fname = 'uniprot_entries.tsv'
    logging.info('Saving into %s.' % fname)
    with open(fname, 'wb') as fh:
        fh.write(full_table.encode('utf-8'))


if __name__ == '__main__':
    download_phosphositeplus()
    #download_uniprot_mappings()
