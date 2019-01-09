from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from nose.plugins.attrib import attr
from protmapper.phosphosite_client import map_to_human_site, sites_only


def test_map_mouse_to_human():
    mouse_up_id = 'Q61337'
    pos = map_to_human_site(mouse_up_id, 'S', '112')
    assert pos == '75'


def test_isoform_mapping_from_human():
    up_id = 'P29353'
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


def test_isoform_mapping_from_mouse():
    up_id = 'P29353'
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


def test_mapping_from_mouse_isoform():
    up_id = 'Q8CI51-3'
    pos = map_to_human_site(up_id, 'S', '105')
    assert pos == '214'


def test_data_as_triples():
    sites = sites_only()
    assert ('Q8CI51-3', 'S', '105') in sites
    assert ('Q8CI51', 'T', '80') in sites
    # T80 exists in mouse isoform 3 as well, but we test the assumption that
    # Phosphosite doesn't include sites for isoforms with the same position
    # as in the reference isoform
    assert ('Q8CI51-3', 'T', '80') in sites
