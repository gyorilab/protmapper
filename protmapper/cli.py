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


if __name__ == '__main__':
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
    args = parser.parse_args()

    pm = ProtMapper()
    sites = process_input(args.input)
    mapped_sites = pm.map_sitelist_to_human_ref(sites)
    dump_output(args.output, mapped_sites)
