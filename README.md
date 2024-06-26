# Protmapper
The Protmapper maps protein sites to the human reference
sequence based on UniProt, PhosphoSitePlus, and manual curation.


## Installation and usage

### Python package
The Protmapper is a Python package that is available on PyPI and can be
installed as:

```bash
pip install protmapper
```

### Docker container
Protmapper can be run as a local service via a Docker container exposing a
REST API as:

```bash
docker run -d -p 8008:8008 gyorilab/protmapper:latest
```

Example: once the container is running, you can send requests to the REST API

```bash
curl -X POST -H "Content-Type: application/json" -d '{"site_list": [["P28482", "uniprot", "T", "184"]]}' http://localhost:8008/map_sitelist_to_human_ref
```

which is equivalent to the following Python code using the `requests` package

```python
import requests
url = 'http://localhost:8008/map_sitelist_to_human_ref'
data = {'site_list': [['P28482', 'uniprot', 'T', '184']]}
response = requests.post(url, json=data)
print(response.json())
```

### Command line interface
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
                        file. Each line represents a mapped site produced by
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

Example: the sample file [cli_input.csv](https://raw.githubusercontent.com/gyorilab/protmapper/master/protmapper/tests/cli_input.csv)
has the following content

```csv
MAPK1,hgnc,T,183
MAPK1,hgnc,T,184
MAPK1,hgnc,T,185
MAPK1,hgnc,T,186
```

By running the following command

```bash
protmapper cli_input.csv output.csv
```

we get `output.csv` which has the following content

```csv
up_id,error_code,valid,orig_res,orig_pos,mapped_id,mapped_res,mapped_pos,description,gene_name
P28482,,False,T,183,P28482,T,185,INFERRED_MOUSE_SITE,MAPK1
P28482,,False,T,184,P28482,T,185,INFERRED_METHIONINE_CLEAVAGE,MAPK1
P28482,,True,T,185,,,,VALID,MAPK1
P28482,,False,T,186,,,,NO_MAPPING_FOUND,MAPK1
```



## Documentation
For a detailed documentation of the Protmapper, visit http://protmapper.readthedocs.io

## Funding
The development of Protmapper is funded under the DARPA grants W911NF018-1-0124
and HR00112220036.

## Citation

```bibtex
@article{bachman2022protmapper,
  author = {Bachman, John A and Sorger, Peter K and Gyori, Benjamin M},
  doi = {10.1101/822668},
  journal = {bioRxiv},
  publisher = {Cold Spring Harbor Laboratory},
  title = {{Assembling a corpus of phosphoproteomic annotations using ProtMapper to normalize site information from databases and text mining}},
  url = {https://www.biorxiv.org/content/10.1101/822668v4},
  year = {2022}
}
```
