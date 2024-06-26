from protmapper.phosphosite_client import map_to_human_site, sites_only, \
                                          PspMapping

def test_map_mouse_to_human():
    mouse_up_id = 'Q61337'
    psp = map_to_human_site(mouse_up_id, 'S', '112')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q92934' # Human ref seq
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '75'
    assert psp.motif == 'EIRSRHSSYPAGTED'
    assert psp.respos == 7


def test_isoform_mapping_from_human():
    up_id = 'P29353'
    psp = map_to_human_site(up_id, 'Y', '239')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P29353' # Human ref seq
    assert psp.mapped_res == 'Y'
    assert psp.mapped_pos == '349'
    assert psp.motif == 'EEPPDHQYYNDFPGK'
    assert psp.respos == 7


def test_mapping_from_human_isoform():
    up_id = 'P29353-2'
    psp = map_to_human_site(up_id, 'Y', '239')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P29353' # Human ref seq
    assert psp.mapped_res == 'Y'
    assert psp.mapped_pos == '349'
    assert psp.motif == 'EEPPDHQYYNDFPGK'
    assert psp.respos == 7


def test_isoform_mapping_from_mouse():
    up_id = 'P98083' # Mouse SHC1
    psp = map_to_human_site(up_id, 'Y', '239')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P29353' # Human ref seq
    assert psp.mapped_res == 'Y'
    assert psp.mapped_pos == '349'
    assert psp.motif == 'EEPPDHQYYNDFPGK'
    assert psp.respos == 7


def test_mapping_from_human_ref():
    up_id = 'P29353' # SHC1
    psp = map_to_human_site(up_id, 'Y', '349')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P29353' # Human ref seq
    assert psp.mapped_res == 'Y'
    assert psp.mapped_pos == '349'
    assert psp.motif == 'EEPPDHQYYNDFPGK'
    assert psp.respos == 7


def test_mapping_from_human_ref_iso_id():
    up_id = 'P29353-1' # SHC1
    psp = map_to_human_site(up_id, 'Y', '349')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P29353' # Human ref seq
    assert psp.mapped_res == 'Y'
    assert psp.mapped_pos == '349'
    assert psp.motif == 'EEPPDHQYYNDFPGK'
    assert psp.respos == 7


def test_mapping_from_mouse_isoform():
    up_id = 'Q8CI51-3'
    psp = map_to_human_site(up_id, 'S', '105')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q96HC4' # Human ref seq
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '214'
    assert psp.motif == 'PTVTSVCSETSQELA'
    assert psp.respos == 7


def test_no_site_in_human_ref():
    psp = map_to_human_site('Q01105', 'S', '9')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q01105-2'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '9'
    assert psp.motif == 'SAPAAKVSKKELNSN'
    assert psp.respos == 7


"""
def test_wrong_residue():
    # SL6A3 T53 -> S53
    psp = map_to_human_site('Q01959', 'T', '53')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q01959'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '53'
    assert psp.motif == 'TLTNPRQSPVEAQDR'
    assert psp.respos == 7
"""


def test_smpd1_s508():
    # The site is invalid, but PSP doesn't know that
    psp = map_to_human_site('P17405', 'S', '508')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'P17405'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '508'
    assert psp.motif == 'DGNYSGSSHVVLDHE'
    assert psp.respos == 7


def test_set_s9():
    psp = map_to_human_site('Q01105', 'S', '9')
    assert isinstance(psp, PspMapping)
    assert psp.mapped_id == 'Q01105-2'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '9'


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


def test_explicit_ref_isoforms():
    psp = map_to_human_site('Q9Y2K2', 'S', '551')
    assert psp.mapped_id == 'Q9Y2K2'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '493'

    psp = map_to_human_site('Q14155', 'S', '672')
    assert psp.mapped_id == 'Q14155'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '694'

    psp = map_to_human_site('O15027', 'T', '220')
    assert psp.mapped_id == 'O15027'
    assert psp.mapped_res == 'T'
    assert psp.mapped_pos == '415'

    psp = map_to_human_site('Q16555', 'S', '627')
    assert psp.mapped_id == 'Q16555'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '522'


def test_ref_seq_not_found():
    psp = map_to_human_site('P10636', 'S', '202')
    assert psp.mapped_id == 'P10636'
    assert psp.mapped_res == 'S'
    assert psp.mapped_pos == '519'

