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
This endpoint takes a single `site_list` argument which is a list of lists
where each list consists of exactly 4 elements in the following order:
`prot_id`, `prot_ns`, `residue`, and `position`. The response
is a list of MappedSite object represented as JSON.

Example
~~~~~~~

Input:

.. code-block:: javascript

    {"site_list": [
        ["MAPK1","hgnc","T","185"],
        ["MAPK1", "hgnc", "T", "183"]
        ]
    }

Output:

.. code-block:: javascript

    [
     {
      "description": "VALID",
      "error_code": null,
      "gene_name": "MAPK1",
      "mapped_id": null,
      "mapped_pos": null,
      "mapped_res": null,
      "orig_pos": "185",
      "orig_res": "T",
      "up_id": "P28482",
      "valid": true
     },
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
    ]

Optional arguments
------------------
Both endpoints take the following optional boolean arguments which are `true`
by default:

- `do_methionine_offset`
- `do_orthology_mapping`
- `do_isoform_mapping`

"""