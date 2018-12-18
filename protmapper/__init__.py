import logging
from sitemapper.api import SiteMapper, MappedSite
from sitemapper.resources import resource_dir


logging.basicConfig(format=('%(levelname)s: [%(asctime)s] %(name)s'
                            ' - %(message)s'),
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
