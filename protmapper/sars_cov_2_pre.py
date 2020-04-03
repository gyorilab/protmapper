import logging
import argparse
from pathlib import Path
from xml.etree import ElementTree as ET

logger = logging.getLogger(__name__)

SARS_NAME = 'SARS-CoV-2'
UP_NS = '{http://uniprot.org/uniprot}'


def process_entry(entry):
    """Process one of the entry tags in the xml file

    Parameters
    ----------
    entry : xml.etree.ElementTree.Element

    Returns
    -------
    list
        A list of name,UP-ID,organism tuples
    """
    # Initialize list
    name_mapping = []

    # Get the UP ID
    up_id = entry.find(UP_NS + 'accession').text
    logger.info('Processing uniprot id %s' % up_id)

    # NOTE: skip the <name> tag for now
    # Get the <name> tag
    # name_ = entry.find(UP_NS + 'name').text
    # name_mapping.append((name_, up_id, SARS_NAME))

    # Get all names:
    # protein -> recommendedName; alternativeName
    #            recommendedName -> fullName; shortName
    #            alternativeName -> fullName; shortName
    protein = entry.find(UP_NS + 'protein')
    for child_tag in protein:
        if child_tag.tag.lower() in {UP_NS + 'recommendedname',
                                     UP_NS + 'alternativename'}:
            for fullname_tag in child_tag.findall(UP_NS + 'fullName'):
                name_mapping.append((fullname_tag.text, up_id, SARS_NAME))
            for shortname_tag in child_tag.findall(UP_NS + 'shortName'):
                name_mapping.append((shortname_tag.text, up_id, SARS_NAME))

    return name_mapping


def process_xml(fname):
    # Read file into xml.etree.ElementTree
    et = ET.parse(fname)

    # Process xml
    name_mappings = []
    for entry in et.findall(UP_NS + 'entry'):
        name_mappings.extend(process_entry(entry))

    return name_mappings
