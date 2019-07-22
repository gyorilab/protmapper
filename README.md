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

## Documentation
For a detailed documentation of the Protmapper, visit http://protmapper.readthedocs.io

## Funding
The development of protmapper is funded under the DARPA Automated Scientific Discovery Framework project (ARO grant W911NF018-1-0124).
