"""The Protmapper REST API allows interacting with the Protmapper through HTTP
requests. The REST API takes GET or POST request with a JSON payload.

The REST API exposes the following endpoints:

`map_to_human_ref`
------------------

This endpoint takes 4 arguments: `prot_id`, `prot_ns`, `residue`, and
`position` and returns a JSON representation of a MappedSite object.

Example
~~~~~~~

Input:

.. code-block:: javascript

    {"prot_id": "MAPK1",
     "prot_ns": "hgnc",
     "residue": "T",
     "position": "183"}

Output:

.. code-block:: javascript

    {
     "description": "INFERRED_MOUSE_SITE",
     "error_code": null,
     "gene_name": "MAPK1",
     "mapped_id": "P28482",
     "mapped_pos": "185",
     "mapped_res": "T",
     "orig_pos": "183",
     "orig_res": "T",
     "up_id": "P28482",
     "valid": false
    }

`map_sitelist_to_human_ref`
---------------------------
This endpoint takes a list of lists as argument where each list consists of
exactly 4 elements in the following order: `prot_id`, `prot_ns`, `residue`,
and `position`.

Example
~~~~~~~

Inp

"""