from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import os
import pickle
from nose.tools import raises
from sitemapper.api import SiteMapper, _validate_site, MappedSite


@raises(ValueError)
def test_invalid_residue():
    sm = SiteMapper()
    sm.map_to_human_ref('MAPK1', 'HGNC', 'B', '185')


@raises(ValueError)
def test_invalid_position():
    sm = SiteMapper()
    sm.map_to_human_ref('MAPK1', 'HGNC', 'T', 'foo')


@raises(ValueError)
def test_invalid_prot_ns():
    sm = SiteMapper()
    sm.map_to_human_ref('MAPK1', 'hgncsymb', 'T', '185')


def test_validate_site():
    valid_site = _validate_site('T', '185')
    assert isinstance(valid_site, tuple)
    assert valid_site[0] == 'T'
    assert valid_site[1] == '185'

    valid_site = _validate_site('T', 185)
    assert isinstance(valid_site, tuple)
    assert valid_site[0] == 'T'
    assert valid_site[1] == '185'


def test_check_agent_mod_up_id():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '185')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'P28482'
    assert ms.gene_name == 'MAPK1'
    assert ms.valid is True
    assert ms.orig_res == 'T'
    assert ms.orig_pos == '185'
    assert ms.mapped_res == 'T'
    assert ms.mapped_pos == '185'
    assert ms.description == 'VALID'

    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '183')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'P28482'
    assert ms.gene_name == 'MAPK1'
    assert ms.valid is False
    assert ms.orig_res == 'T'
    assert ms.orig_pos == '183'
    assert ms.mapped_res == 'T'
    assert ms.mapped_pos == '185'
    assert ms.description == 'INFERRED_MOUSE_SITE'


def test_check_agent_mod_hgnc():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '185')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'P28482'
    assert ms.gene_name == 'MAPK1'
    assert ms.valid is True
    assert ms.orig_res == 'T'
    assert ms.orig_pos == '185'
    assert ms.mapped_res == 'T'
    assert ms.mapped_pos == '185'
    assert ms.description == 'VALID'

    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '183')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'P28482'
    assert ms.gene_name == 'MAPK1'
    assert ms.valid is False
    assert ms.orig_res == 'T'
    assert ms.orig_pos == '183'
    assert ms.mapped_res == 'T'
    assert ms.mapped_pos == '185'
    assert ms.description == 'INFERRED_MOUSE_SITE'


def test_map_mouse_site():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('THEMIS2', 'hgnc', 'Y', '660')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'Q5TEJ8'
    assert ms.gene_name == 'THEMIS2'
    assert ms.valid is False
    assert ms.orig_res == 'Y'
    assert ms.orig_pos == '660'
    assert ms.mapped_res == 'Y'
    assert ms.mapped_pos == '632'
    assert ms.description == 'INFERRED_MOUSE_SITE'


def test_map_rat_site():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('NPHS1', 'hgnc', 'Y', '1204')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'O60500'
    assert ms.gene_name == 'NPHS1'
    assert ms.valid is False
    assert ms.orig_res == 'Y'
    assert ms.orig_pos == '1204'
    assert ms.mapped_res == 'Y'
    assert ms.mapped_pos == '1193'
    assert ms.description == 'INFERRED_RAT_SITE'


def test_map_methionine_cleavage():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('DAXX', 'hgnc', 'S', '667')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'Q9UER7'
    assert ms.gene_name == 'DAXX'
    assert ms.valid is False
    assert ms.orig_res == 'S'
    assert ms.orig_pos == '667'
    assert ms.mapped_res == 'S'
    assert ms.mapped_pos == '668'
    assert ms.description == 'INFERRED_METHIONINE_CLEAVAGE'

def test_repr_str():
    sm = SiteMapper()
    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '183')
    assert str(ms) == ("MappedSite(up_id='P28482', valid=False, "
                       "orig_res='T', orig_pos='183', mapped_res='T', "
                       "mapped_pos='185', description='INFERRED_MOUSE_SITE', "
                       "gene_name='MAPK1')")

def test_read_cache():
    sm = SiteMapper(use_cache=True, cache_path='test_cache_read.pkl')
    ms = sm.map_to_human_ref('P28482', 'uniprot', 'Q', '32')
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'TEST1'
    assert ms.gene_name == 'TEST2'
    assert ms.valid == 'TEST3'
    assert ms.orig_res == 'TEST4'
    assert ms.orig_pos == 'TEST5'
    assert ms.mapped_res == 'TEST6'
    assert ms.mapped_pos == 'TEST7'
    assert ms.description == 'TEST8'

def test_write_cache():
    cache_path = 'test_cache_write.pkl'
    assert not os.path.isfile(cache_path)
    sm = SiteMapper(use_cache=True, cache_path='test_cache_write.pkl')
    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '183')
    del sm
    assert os.path.isfile(cache_path)
    with open(cache_path, 'rb') as f:
        cache_dict = pickle.load(f)
    assert isinstance(cache_dict, dict)
    assert len(cache_dict) == 1
    ms = cache_dict[('P28482', 'T', '183')]
    assert isinstance(ms, MappedSite)
    assert ms.up_id == 'P28482'
    assert ms.gene_name == 'MAPK1'
    assert ms.valid is False
    assert ms.orig_res == 'T'
    assert ms.orig_pos == '183'
    assert ms.mapped_res == 'T'
    assert ms.mapped_pos == '185'
    assert ms.description == 'INFERRED_MOUSE_SITE'
    os.unlink(cache_path)


def test_mapped_site_equals():
    ms1 = MappedSite(up_id='P28482', valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms2 = MappedSite(up_id='P28482', valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    assert ms1 == ms2
    ms2.gene_name = 'FOO'
    assert ms1 != ms2

def test_mapped_site_hash():
    ms1 = MappedSite(up_id='P28482', valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms2 = MappedSite(up_id='P28482', valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    assert hash(ms1) == hash(ms2)
    ms2.gene_name = 'FOO'
    assert hash(ms1) != hash(ms2)


if __name__ == '__main__':
    test_mapped_site_hash()
