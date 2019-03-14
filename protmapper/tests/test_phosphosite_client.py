from nose.plugins.attrib import attr
from protmapper.phosphosite_client import map_to_human_site, sites_only, \
                                          PspMapping

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


def test_no_site_in_human_ref():
    psp = map_to_human_site('Q01105', 'S', '9')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q01105-2'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '9'
    assert psp.motif == 'SAPAAKVSKKELNSN'
    assert psp.respos == 7


#def test_wrong_residue():
#    # SL6A3 T53 -> S53
#    psp = map_to_human_site('Q01959', 'T', '53')
#    assert isinstance(psp, PspMapping)
#    assert psp.mapped_id == 'Q01959'
#    assert psp.mapped_res == 'S'
#    assert psp.mapped_pos == '53'
#    assert psp.motif == 'TLTNPRQSPVEAQDR'
#    assert psp.respos == 7


def test_smpd1_s508():
    # The site is invalid, but PSP doesn't know that
    psp = map_to_human_site('P17405', 'S', '508')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P17405'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '508'
    assert psp.motif == 'DGNYSGSSHVVLDHE'
    assert psp.respos == 7


def test_h2afx_s139():
    psp = map_to_human_site('P16104', 'S', '139')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P16104'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '139'
    assert psp.motif == 'GKKATQASQEY'
    assert psp.respos == 7


def test_motif_processing():
    # Make sure that site motifs with prepended underscores have the residue
    # position assigned accordingly
    psp = map_to_human_site('P68431', 'T', '3')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P68431'
    assert psp.mapped_res == 'T'
    assert psp.mapped_pos == '3'
    assert psp.motif == 'ARTKQTARKS'
    assert psp.respos == 2


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


if __name__ == '__main__':
    test_motif_processing()
    test_no_site_in_human_ref()
    test_smpd1_s508()
    test_h2afx_s139()
