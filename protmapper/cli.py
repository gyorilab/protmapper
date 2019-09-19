import csv
import argparse
from protmapper.api import *


def process_input(fname):
    sites = []
    with open(fname, 'r') as fh:
        for idx, row in enumerate(csv.reader(fh)):
            if len(row) != 4:
                raise ValueError('Line %d of %s doesn\'t have 4 elements.' %
                                 (idx, fname))
            sites.append(row)
    return sites


def dump_output(fname, mapped_sites):
    rows = [mapped_sites[0].attrs]
    rows += [ms.to_list() for ms in mapped_sites]
    with open(fname, 'w') as fh:
        writer = csv.writer(fh)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(
        description='Run Protmapper on a list of proteins with residues and '
                    'sites provided in a text file.')
    parser.add_argument('input',
        help=('Path to an input file. The input file is a text file in '
              'which each row consists of four comma separated values, '
              'with the first element being a protein ID, the second, '
              'the namespace in which that ID is valid (uniprot or hgnc),'
              'third, an amino acid represented as a single capital letter, '
              'and fourth, a site position on the protein.'))
    parser.add_argument('output',
        help=('Path to the output file to be generated. Each line of the '
              'output file corresponds to a line in the input file. Each line'
              'represents a mapped site produced by Protmapper.'))
    parser.add_argument('--peptide',
        help=('If given, the third element of each row of the input file is a '
              'peptide (amino acid sequence) rather than a single amino acid '
              'residue. In this case, peptide-oriented mappings are '
              'applied. In this mode the following boolean arguments are '
              'ignored.'), action='store_true')
    parser.add_argument('--no_methionine_offset', help=(
        'If given, will not check for off-by-one errors in site position ('
        'possibly) attributable to site numbering from mature proteins after '
        'cleavage of the initial methionine.'), action='store_true')
    parser.add_argument('--no_orthology_mapping', help=(
        'If given, will not check sequence positions for known modification '
        'sites in mouse or rat sequences (based on PhosphoSitePlus data).'),
        action='store_true')
    parser.add_argument('--no_isoform_mapping', help=(
        'If given, will not check sequence positions for known modifications '
        'in other human isoforms of the protein (based on PhosphoSitePlus '
        'data).'), action='store_true')
    args = parser.parse_args()
    # Separate call to make function testable
    run_main(args)


def run_main(args):
    mapping_kwargs = {
        'do_methionine_offset': False if args.no_methionine_offset else True,
        'do_orthology_mapping': False if args.no_orthology_mapping else True,
        'do_isoform_mapping': False if args.no_isoform_mapping else True
    }

    pm = ProtMapper()
    sites = process_input(args.input)
    if args.peptide:
        # We have to make the positions ints here
        sites = [tuple(s[:3] + [int(s[3])]) for s in sites]
        mapped_sites = [pm.map_peptide_to_human_ref(*site) for site in sites]
    else:
        mapped_sites = pm.map_sitelist_to_human_ref(sites, **mapping_kwargs)
    dump_output(args.output, mapped_sites)


if __name__ == '__main__':
    main()
