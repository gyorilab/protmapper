from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
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

