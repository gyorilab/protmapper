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


def test_mapping_from_human_isoform():
    up_id = 'P29353-2'
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


def test_isoform_mapping_from_mouse():
    up_id = 'P29353' # SHC1
    pos = map_to_human_site(up_id, 'Y', '239')
    assert pos == '349'


def test_mapping_from_human_ref():
    up_id = 'P29353' # SHC1
    pos = map_to_human_site(up_id, 'Y', '349')
    assert pos == '349'


def test_mapping_from_human_ref_iso_id():
    up_id = 'P29353-1' # SHC1
    pos = map_to_human_site(up_id, 'Y', '349')
    assert pos == '349'


def test_mapping_from_mouse_isoform():
    up_id = 'Q8CI51-3'
    pos = map_to_human_site(up_id, 'S', '105')
    assert pos == '214'


def test_sites_only():
    sites = sites_only()
    # These checks make sure that the sites are constructed from the data
    # dictionaries as expected
    assert ('Q8CI51-3', 'S', '105') in sites
    assert ('Q8CI51', 'T', '80') in sites
    assert ('Q8CI51-3', 'T', '80') not in sites
    assert ('P29353', 'Y', '427') in sites
    assert ('P29353', 'Y', '317') in sites
    assert ('P29353-2', 'Y', '317') in sites
    assert ('P29353-2', 'Y', '427') not in sites
    # Check the -1 isoforms are also included
    assert ('P28661-1', 'S', '28') in sites

