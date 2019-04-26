__version__ = '0.0.7'
import logging

logging.basicConfig(format=('%(levelname)s: [%(asctime)s] %(name)s'
                            ' - %(message)s'),
                    level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')

logging.getLogger('requests').setLevel(logging.ERROR)
logging.getLogger('urllib3').setLevel(logging.ERROR)
logging.getLogger('rdflib').setLevel(logging.ERROR)
logging.getLogger('boto3').setLevel(logging.CRITICAL)
logging.getLogger('botocore').setLevel(logging.CRITICAL)

logger = logging.getLogger('protmapper')

from protmapper.api import ProtMapper, MappedSite
from protmapper.resources import resource_dir
