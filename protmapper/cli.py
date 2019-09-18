import csv
import argparse
from protmapper.api import *


def process_input(fname):
    sites = []
    with open(fname, 'r') as fh:
        for idx, row in enumerate(csv.reader(fh)):
            if len(row) != 3:
                raise ValueError('Line %d doesn\'t have 3 elements.')
            sites.append(row)
    return sites


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run Protmapper on a list of proteins with residues and '
                    'sites provided in a text file.')
    parser.add_argument('--input',
        help=('Path to an input file. The input file is a text file in '
              'which each row consists of four comma separated values, '
              'with the first element being a protein ID, the second, '
              'the namespace in which that ID is valid (e.g., UP, HGNC),'
              'third, an amino acid represented as a single capital letter, '
              'and fourth, a site position on the protein.'), required=True)
    parser.add_argument('--output',
        help=('Path to the output file to be generated. Each line of the '
              'output file corresponds to a line in the input file. Each line'
              'represents a mapped site produced by Protmapper.'),
              required=True)
    args = parser.parse_args()

    pm = ProtMapper()
    sites = process_input(args.input)
    mapped_sites = pm.map_sitelist_to_human_ref(sites)