import os
from collections import Counter
from os.path import join, abspath, dirname, isfile
import pickle
from nose.tools import raises
from protmapper.api import ProtMapper, _validate_site, MappedSite, \
                           InvalidSiteException


@raises(InvalidSiteException)
def test_validate_invalid_residue1():
    _validate_site('B', '185')


@raises(InvalidSiteException)
def test_validate_invalid_residue2():
    _validate_site('T', 'foo')


@raises(InvalidSiteException)
def test_validate_invalid_residue3():
    _validate_site('T', '12.5')


def test_validate_site():
    valid_site = _validate_site('T', '185')
    assert isinstance(valid_site, tuple)
    assert valid_site[0] == 'T'
    assert valid_site[1] == '185'

    valid_site = _validate_site('T', 185)
    assert isinstance(valid_site, tuple)
    assert valid_site[0] == 'T'
    assert valid_site[1] == '185'


def test_mapped_site_equals():
    ms1 = MappedSite(up_id='P28482', error_code=None, valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms2 = MappedSite(up_id='P28482', valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    assert ms1 == ms2
    ms2.gene_name = 'FOO'
    assert ms1 != ms2
    assert ms1.not_invalid() is False
    assert ms1.has_mapping() is True
    assert ms2.not_invalid() is False
    assert ms2.has_mapping() is True


def test_mapped_site_hash():
    ms1 = MappedSite(up_id='P28482', error_code=None, valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms2 = MappedSite(up_id='P28482', error_code=None, valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    assert hash(ms1) == hash(ms2)
    ms2.gene_name = 'FOO'
    assert hash(ms1) != hash(ms2)


def test_mapped_site_set_ctr():
    """Check if two identical sites in different objects are handled as equal
    in sets and counters."""
    ms1 = MappedSite(up_id='P28482', error_code=None, valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms2 = MappedSite(up_id='P28482', error_code=None, valid=False, orig_res='T',
                     orig_pos='183', mapped_res='T', mapped_pos='185',
                     description='INFERRED_MOUSE_SITE', gene_name='MAPK1')
    ms_list = [ms1, ms2]
    assert len(set(ms_list)) == 1
    ctr = Counter(ms_list)
    assert len(ctr) == 1
    assert list(ctr.values())[0] == 2


def test_map_invalid_residue():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('MAPK1', 'hgnc', 'B', '185')
    assert ms == MappedSite('P28482', None, 'B', '185',
                            error_code='INVALID_SITE',
                            description='Residue B not a valid amino acid')
    assert ms.not_invalid() is True
    assert ms.has_mapping() is False


def test_map_invalid_position1():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('MAPK1', 'hgnc', 'T', 'foo')
    assert ms == MappedSite(
            'P28482', None, 'T', 'foo', error_code='INVALID_SITE',
            description='Position foo not a valid sequence position.')


def test_map_invalid_position2():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('MAPK1', 'hgnc', 'T', '12.5')
    assert ms == MappedSite(
            'P28482', None, 'T', '12.5', error_code='INVALID_SITE',
            description='Position 12.5 not a valid sequence position.')


@raises(ValueError)
def test_invalid_prot_ns():
    sm = ProtMapper()
    sm.map_to_human_ref('MAPK1', 'hgncsymb', 'T', '185')


def test_check_agent_mod_up_id():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '185')
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=True,
                            orig_res='T', orig_pos='185', mapped_res=None,
                            mapped_pos=None, description='VALID',
                            gene_name='MAPK1')
    assert ms.not_invalid() is True
    assert ms.has_mapping() is False

    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '183')
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=False,
                            orig_res='T', orig_pos='183', mapped_res='T',
                            mapped_pos='185', description='INFERRED_MOUSE_SITE',
                            gene_name='MAPK1')
    assert ms.not_invalid() is False
    assert ms.has_mapping() is True


def test_check_agent_mod_hgnc():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '185')
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=True,
                            orig_res='T', orig_pos='185', mapped_res=None,
                            mapped_pos=None, description='VALID',
                            gene_name='MAPK1')

    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '183')
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=False,
                            orig_res='T', orig_pos='183', mapped_res='T',
                            mapped_pos='185', description='INFERRED_MOUSE_SITE',
                            gene_name='MAPK1')


def test_map_mouse_site():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('THEMIS2', 'hgnc', 'Y', '660')
    assert ms == MappedSite(up_id='Q5TEJ8', error_code=None, valid=False,
                            orig_res='Y', orig_pos='660', mapped_res='Y',
                            mapped_pos='632', description='INFERRED_MOUSE_SITE',
                            gene_name='THEMIS2')


def test_map_rat_site():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('NPHS1', 'hgnc', 'Y', '1204')
    assert ms == MappedSite(up_id='O60500', error_code=None, valid=False,
                            orig_res='Y', orig_pos='1204', mapped_res='Y',
                            mapped_pos='1193', description='INFERRED_RAT_SITE',
                            gene_name='NPHS1')


def test_map_methionine_cleavage():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('DAXX', 'hgnc', 'S', '667')
    assert ms == MappedSite(up_id='Q9UER7', error_code=None, valid=False,
                            orig_res='S', orig_pos='667', mapped_res='S',
                            mapped_pos='668',
                            description='INFERRED_METHIONINE_CLEAVAGE',
                            gene_name='DAXX')


def test_map_from_sitemap():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('P15056', 'uniprot', 'S', '753')
    assert ms == MappedSite(up_id='P15056', valid=False, orig_res='S',
                            orig_pos='753', mapped_res='T', mapped_pos='753',
                            description='wrong residue', gene_name='BRAF')


def test_map_invalid_from_sitemap():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('P06400', 'uniprot', 'D', '1')
    assert ms == MappedSite(up_id='P06400', valid=False, orig_res='D',
                            orig_pos='1', mapped_res=None, mapped_pos=None,
                            description='probable reading error: cyclin D1',
                            gene_name='RB1')


def test_h2afx_s139():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('P16104', 'uniprot', 'S', '139')
    assert ms == MappedSite(up_id='P16104', valid=False, orig_res='S',
                            orig_pos='139', mapped_res='S', mapped_pos='140',
                            description='SEQ_MISMATCH_PSP_UP',
                            gene_name='H2AFX')


def test_sl6a3_t53():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('Q01959', 'uniprot', 'T', '53')
    assert ms == MappedSite(up_id='Q01959', valid=False, orig_res='T',
                            orig_pos='53', mapped_res='S', mapped_pos='53',
                            description='INFERRED_WRONG_RESIDUE',
                            gene_name='Q01959')


def test_smpd1_s508():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('P17405', 'uniprot', 'S', '508')
    assert ms == MappedSite(up_id='P17405', valid=False, orig_res='S',
                            orig_pos='508', mapped_res='S', mapped_pos='510',
                            description='INFERRED_ALTERNATIVE_ISOFORM',
                            gene_name='SMPD1')


def test_set_s9():
    pm = ProtMapper()
    ms = pm.map_to_human_ref('Q01105', 'uniprot', 'S', '9')
    assert ms == MappedSite(up_id='Q01105', valid=False, orig_res='S',
                            orig_pos='9', mapped_res='S', mapped_pos='9',
                            description='INFERRED_ALTERNATIVE_ISOFORM',
                            gene_name='SET')


def test_repr_str():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('MAPK1', 'hgnc', 'T', '183')
    assert str(ms) == ("MappedSite(up_id='P28482', error_code=None, "
                       "valid=False, orig_res='T', orig_pos='183', "
                       "mapped_res='T', mapped_pos='185', "
                       "description='INFERRED_MOUSE_SITE', "
                       "gene_name='MAPK1')")


def test_write_cache():
    cache_path = join(dirname(abspath(__file__)), 'test_cache_write.pkl')
    assert not isfile(cache_path)
    sm = ProtMapper(use_cache=True, cache_path=cache_path)
    ms = sm.map_to_human_ref('P28482', 'uniprot', 'T', '183')
    del sm
    assert isfile(cache_path)
    with open(cache_path, 'rb') as f:
        cache_dict = pickle.load(f)
    assert isinstance(cache_dict, dict)
    assert len(cache_dict) == 1
    ms = cache_dict[('P28482', 'T', '183')]
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=False,
                            orig_res='T', orig_pos='183', mapped_res='T',
                            mapped_pos='185',
                            description='INFERRED_MOUSE_SITE',
                            gene_name='MAPK1')
    os.unlink(cache_path)


def test_mapped_site_to_json():
    kwargs = dict(up_id='P28482', error_code=None, valid=True, orig_res='T',
                orig_pos='185', mapped_res='T', mapped_pos='185',
                description='VALID', gene_name='MAPK1')
    ms = MappedSite(**kwargs)
    assert ms.to_json() == kwargs


def test_invalid_uniprot_http_error():
    """Trigger an HTTP Error by passing in an invalid Uniprot ID. The error
    should be indicated in the returned MappedSite."""
    sm = ProtMapper()
    ms = sm.map_to_human_ref('ASDF', 'uniprot', 'Q', '999')
    assert ms.error_code == 'UNIPROT_HTTP_NOT_FOUND'
    assert ms == MappedSite(up_id='ASDF', error_code='UNIPROT_HTTP_NOT_FOUND',
                   valid=None, orig_res='Q', orig_pos='999',
                   mapped_res=None, mapped_pos=None, description=None,
                   gene_name=None)


def test_invalid_gene_name_error():
    sm = ProtMapper()
    ms = sm.map_to_human_ref('ASDF', 'hgnc', 'Q', '999')
    assert ms.error_code == 'NO_UNIPROT_ID'
    assert ms == MappedSite(up_id=None, error_code='NO_UNIPROT_ID',
                   valid=None, orig_res='Q', orig_pos='999',
                   mapped_res=None, mapped_pos=None, description=None,
                   gene_name='ASDF')


def test_motif_from_position():
    sm = ProtMapper()
    # Try with string position
    (motif, respos) = sm.motif_from_position('O95218-2', '120', window=7)
    assert motif == 'EYIEREESDGEYDEF'
    assert respos == 8
    # Also try with int
    motif, respos = sm.motif_from_position('O95218-2', '120')
    assert motif == 'EYIEREESDGEYDEF'
    assert respos == 8
    # Change the window size
    motif, respos = sm.motif_from_position('O95218-2', '120', window=5)
    assert motif == 'IEREESDGEYD'
    assert respos == 6
    # Try with a residue that is close to the end of the protein
    motif, respos = sm.motif_from_position('Q04637-1', '1596', window=7)
    assert motif == 'LREAEEESDHN'
    assert respos == 8
    # Try with a residue that is close to the start of the protein
    motif, respos = sm.motif_from_position('Q04637-1', '3', window=7)
    assert motif == 'MNKAPQSTGP'
    assert respos == 3
    (motif, respos) = sm.motif_from_position('P28482', '187')
    assert motif == 'HTGFLTEYVATRWYR'
    assert respos == 8


def test_map_peptide():
    sm = ProtMapper()
    # Check the localization of a peptide to a target protein
    peptide = 'MNTPSQPRQHFY'
    pos = 5
    respos = sm.map_peptide('Q04637-1', peptide, pos)
    assert respos == 45
    # Check that the motif is not found in an unrelated sequence
    respos = sm.map_peptide('O95218-2', peptide, pos)
    assert respos is None


def test_map_peptide_to_human_ref():
    sm = ProtMapper()
    peptide = 'MNTPSQPRQHFY'
    pos = 5
    ms = sm.map_peptide_to_human_ref('Q04637', 'uniprot', peptide, pos)
    assert isinstance(ms, MappedSite)
    assert ms == MappedSite(up_id='Q04637', error_code=None, valid=True,
                            orig_res=None, orig_pos=None, mapped_res='S',
                            mapped_pos='45',
                            description=None,
                            gene_name='EIF4G1')
    # The same as above except with the gene name instead
    ms2 = sm.map_peptide_to_human_ref('EIF4G1', 'hgnc', peptide, pos)
    assert ms2 == ms


def test_map_peptide_to_human_ref2():
    pm = ProtMapper()
    up_id = 'P07942'
    peptide = 'GDNLLDSRMEIRE'
    sitepos = 7
    ms = pm.map_peptide_to_human_ref(up_id, 'uniprot', peptide, sitepos)
    assert ms == MappedSite(up_id='P07942', error_code=None, valid=True,
                            orig_res=None, orig_pos=None, mapped_res='S',
                            mapped_pos='250', description=None,
                            gene_name='LAMB1')
    assert ms


def test_peptide_round_trip():
    pm = ProtMapper()
    pos = '187'
    motif, site_pos = pm.motif_from_position('P28482', '187')
    ms = pm.map_peptide_to_human_ref('MAPK1', 'hgnc', motif, site_pos)
    assert isinstance(ms, MappedSite)
    assert ms == MappedSite(up_id='P28482', error_code=None, valid=True,
                            orig_res=None, orig_pos=None, mapped_res='Y',
                            mapped_pos='187',
                            description=None,
                            gene_name='MAPK1')
