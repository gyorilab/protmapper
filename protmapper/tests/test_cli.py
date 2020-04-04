import os
import csv
from protmapper.cli import run_main


here = os.path.abspath(os.path.dirname(__file__))


class ArgparseMock(object):
    def __init__(self, input, output, peptide=None, no_methionine_offset=None,
                 no_orthology_mapping=None, no_isoform_mapping=None):
        self.input = input
        self.output = output
        self.peptide = peptide
        self.no_methionine_offset = no_methionine_offset
        self.no_orthology_mapping = no_orthology_mapping
        self.no_isoform_mapping = no_isoform_mapping


def test_default():
    input = os.path.join(here, 'cli_input.csv')
    output = os.path.join(here, 'test_output.csv')
    args = ArgparseMock(input, output)
    run_main(args)
    with open(output, 'r') as fh:
        rows = [r for r in csv.reader(fh)]
    assert rows[1][-2] == 'INFERRED_MOUSE_SITE', rows[1]
    assert rows[2][-2] == 'INFERRED_METHIONINE_CLEAVAGE', rows[2]
    assert rows[3][-2] == 'VALID', rows[3]
    assert rows[4][-2] == 'NO_MAPPING_FOUND', rows[4]


def test_options():
    input = os.path.join(here, 'cli_input.csv')
    output = 'test_cli_output.csv'
    args = ArgparseMock(input, output, no_methionine_offset=True)
    run_main(args)
    with open(output, 'r') as fh:
        rows = [r for r in csv.reader(fh)]
    assert rows[1][-2] == 'INFERRED_MOUSE_SITE', rows[1]
    assert rows[2][-2] == 'NO_MAPPING_FOUND', rows[2]
    assert rows[3][-2] == 'VALID', rows[3]
    assert rows[4][-2] == 'NO_MAPPING_FOUND', rows[4]


def test_peptide():
    input = os.path.join(here, 'cli_input_peptide.csv')
    output = 'test_cli_output.csv'
    args = ArgparseMock(input, output, peptide=True)
    run_main(args)
    with open(output, 'r') as fh:
        rows = [r for r in csv.reader(fh)]
    print(rows)
    assert rows[1][2] == 'True', rows
    assert rows[2][2] == 'False'
