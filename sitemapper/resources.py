from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import logging
import requests


logger = logging.getLogger(__name__)


# If the sitemapper resource directory does not exist, try to create it
home_dir = os.path.expanduser('~')
resource_dir = os.path.join(home_dir, '.sitemapper')


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


if __name__ == '__main__':
    download_phosphositeplus()
