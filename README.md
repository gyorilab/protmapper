# Protmapper
The Protmapper maps references to protein sites to the human reference
sequence based on UniProt, PhosphoSitePlus, and manual curation.


## Installation

### Python package
The Protmapper is a Python package that is available on PyPI and can be
installed as:

```
pip install protmapper
```

### Docker container
Alternatively, the Protmapper Docker container can be run to expose it as
a REST API as:

```
docker run -d -p 8008:8008 labsyspharm/protmapper:latest
```

## Command line interface
In addition to supporting usage via a Python API and a REST service,
Protmapper also provides a command line interface that can be used as follows.

```bash
Run Protmapper on a list of proteins with residues and sites provided in a
text file.

positional arguments:
  input                 Path to an input file. The input file is a text file
                        in which each row consists of four comma separated
                        values, with the first element being a protein ID, the
                        second, the namespace in which that ID is valid
                        (uniprot or hgnc),third, an amino acid represented as
                        a single capital letter, and fourth, a site position
                        on the protein.
  output                Path to the output file to be generated. Each line of
                        the output file corresponds to a line in the input
                        file. Each linerepresents a mapped site produced by
                        Protmapper.

optional arguments:
  -h, --help            show this help message and exit
  --peptide             If given, the third element of each row of the input
                        file is a peptide (amino acid sequence) rather than a
                        single amino acid residue. In this case, peptide-
                        oriented mappings are applied. In this mode the
                        following boolean arguments are ignored.
  --no_methionine_offset
                        If given, will not check for off-by-one errors in site
                        position (possibly) attributable to site numbering
                        from mature proteins after cleavage of the initial
                        methionine.
  --no_orthology_mapping
                        If given, will not check sequence positions for known
                        modification sites in mouse or rat sequences (based on
                        PhosphoSitePlus data).
  --no_isoform_mapping  If given, will not check sequence positions for known
                        modifications in other human isoforms of the protein
                        (based on PhosphoSitePlus data).

```

## Documentation
For a detailed documentation of the Protmapper, visit http://protmapper.readthedocs.io

## Funding
The development of protmapper is funded under the DARPA Automated Scientific Discovery Framework project (ARO grant W911NF018-1-0124).
