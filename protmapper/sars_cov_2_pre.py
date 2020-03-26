import wget
import argparse
from bs4 import BeautifulSoup
from pathlib import Path

pre_release_url = \
    'ftp://ftp.uniprot.org/pub/databases/uniprot/pre_release/coronavirus.xml'
SARS_NAME = 'SARS-CoV-2'


def _ftp_download(path='.', url=pre_release_url):
    outdir = Path(path)
    wget.download(url, out=outdir.as_posix())


def process_entry(entry):
    """Process one of the entry tags in the xml file

    Parameters
    ----------
    entry : bs4.element.Tag

    Returns
    -------
    list
        A list of name,UP-ID,organism tuples
    """
    # Initialize list
    name_mapping = []

    # Get the UP ID
    acc_tag = entry.find('accession')
    up_id = acc_tag.text

    # Get the tag after <accession>, the <name> tag
    name_tag = acc_tag.findNextSibling()
    name_mapping.append((name_tag.text, up_id, SARS_NAME))

    # Get all names:
    # protein -> recommendedname; alternativename
    #            recommendedname -> fullname; shortname
    #            alternativename -> fullname; shortname
    protein = entry.protein
    for tag in protein.children:
        if tag.name in {'recommendedname', 'alternativename'}:
            for fullname_tag in tag.findAll('fullname'):
                name_mapping.append((fullname_tag.text, up_id, SARS_NAME))
            for shortname_tag in tag.findAll('shortname'):
                name_mapping.append((shortname_tag.text, up_id, SARS_NAME))

    return name_mapping


def process_xml(fname):
    # Read file into bs4
    with open(fname, 'r') as xmlf:
        soup = BeautifulSoup(xmlf, 'lxml')

    # Process xml
    name_mappings = []
    for entry in soup.findAll('entry'):
        name_mappings.extend(process_entry(entry))

    return name_mappings


def main(tsv_outfile, ftp_path=None):
    # Download file
    fname = SARS_NAME + '_prerelease.xml'
    _ftp_download(fname)

    # Get name mappings
    sars_mappings = process_xml(fname)

    # Write to tsv
    tsv_outfile = tsv_outfile if tsv_outfile.endswith('.tsv') else \
        tsv_outfile.split('.')[0] + '.tsv'
    with open(tsv_outfile, 'w') as tsvf:
        for mapping in sars_mappings:
            tsvf.write('%s\n' % '\t'.join(mapping))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv-out', required=True,
                        help='The file path for the tsv file of the output')
    # parser.add_argument(
    #     '--download-path',
    #     help='The path to where to download the xml resource file'
    # )

    args = parser.parse_args()

    main(args.tsv_out)

