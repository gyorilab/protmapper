from protmapper import uniprot_client
from protmapper.resources import _process_feature
from nose.plugins.attrib import attr


@attr('webservice')
def test_query_protein_exists():
    g = uniprot_client.query_protein('P00533')
    assert g is not None


@attr('webservice')
def test_query_protein_nonexist():
    g = uniprot_client.query_protein('XXXX')
    assert g is None


@attr('webservice')
def test_query_protein_deprecated():
    g = uniprot_client.query_protein('Q8NHX1')
    assert g is not None
    gene_name = uniprot_client.get_gene_name('Q8NHX1')
    assert gene_name == 'MAPK3'
    gene_name = uniprot_client.get_gene_name('Q8NHX1', web_fallback=False)
    assert gene_name == 'MAPK3'


@attr('webservice')
def test_get_family_members():
    members = uniprot_client.get_family_members('RAF')
    assert 'ARAF' in members
    assert 'BRAF' in members
    assert 'RAF1' in members


def test_get_gene_name_human():
    gene_name = uniprot_client.get_gene_name('P00533')
    assert gene_name == 'EGFR'


'''
The below test can be used as a template
to test for entries missing from the resource
file that are available in the web service. It is only relevant if the
resource file is not up to date.

def test_get_gene_name_no_table_entry():
    gene_name = uniprot_client.get_gene_name('P01814', web_fallback=True)
    assert gene_name == 'IGHV2-70'
    gene_name = uniprot_client.get_gene_name('P01814', web_fallback=False)
    assert gene_name is None
'''


def test_get_synonyms():
    upid = 'Q02750'  # This is MAP2K1
    gene_synonyms = uniprot_client.get_gene_synonyms(upid)
    assert gene_synonyms, gene_synonyms
    assert 'MEK1' in gene_synonyms
    protein_synonyms = uniprot_client.get_protein_synonyms(upid)
    assert protein_synonyms, protein_synonyms
    assert 'MKK1' in protein_synonyms
    all_synonyms = uniprot_client.get_synonyms(upid)
    assert set(gene_synonyms + protein_synonyms) == set(all_synonyms)
    assert 'MAP2K1' in all_synonyms


def test_get_gene_name_nonhuman():
    gene_name = uniprot_client.get_gene_name('P31938')
    assert gene_name == 'Map2k1'


def test_get_gene_name_unreviewed():
    gene_name = uniprot_client.get_gene_name('X6RK18', web_fallback=False)
    assert gene_name == 'EXO5'


@attr('webservice')
def test_get_gene_name_no_gene_name():
    gene_name = uniprot_client.get_gene_name('P04434', web_fallback=False)
    assert gene_name is None
    gene_name = uniprot_client.get_gene_name('P04434', web_fallback=True)
    assert gene_name is None


def test_get_gene_name_multiple_gene_names():
    gene_name = uniprot_client.get_gene_name('Q5VWM5')
    assert gene_name == 'PRAMEF9'


def test_is_human():
    assert uniprot_client.is_human('P00533')


def test_not_is_human():
    assert not uniprot_client.is_human('P31938')


def test_noentry_is_human():
    assert not uniprot_client.is_human('XXXX')


@attr('webservice')
def test_get_sequence():
    seq = uniprot_client.get_sequence('P00533')
    assert len(seq) > 1000


@attr('webservice')
def test_get_modifications():
    mods = uniprot_client.get_modifications('P27361')
    assert ('Phosphothreonine', 202) in mods
    assert ('Phosphotyrosine', 204) in mods


@attr('webservice')
def test_verify_location():
    assert uniprot_client.verify_location('P27361', 'T', 202)
    assert not uniprot_client.verify_location('P27361', 'S', 202)
    assert not uniprot_client.verify_location('P27361', 'T', -1)
    assert not uniprot_client.verify_location('P27361', 'T', 10000)


def test_get_mnemonic():
    mnemonic = uniprot_client.get_mnemonic('Q02750')
    assert mnemonic == 'MP2K1_HUMAN'


def test_is_secondary_primary():
    assert not uniprot_client.is_secondary('Q02750')


def test_is_secondary_secondary():
    assert uniprot_client.is_secondary('Q96J62')


def test_get_primary_id_primary():
    assert uniprot_client.get_primary_id('Q02750') == 'Q02750'


def test_get_primary_id_secondary_hashuman():
    assert uniprot_client.get_primary_id('Q96J62') == 'P61978'


def test_get_primary_id_secondary_nohuman():
    assert uniprot_client.get_primary_id('P31848') in \
        ['P0A5M5', 'P9WIU6', 'P9WIU7']


def test_mouse_from_up():
    assert uniprot_client.get_mgi_id('P28028') == '88190'


def test_up_from_mouse():
    assert uniprot_client.get_id_from_mgi('88190') == 'P28028'


def test_rat_from_up():
    assert uniprot_client.get_rgd_id('O08773') == '620003'


def test_up_from_rat():
    assert uniprot_client.get_id_from_rgd('620003') == 'O08773'


def test_mouse_from_human():
    assert uniprot_client.get_mouse_id('P15056') == 'P28028'


def test_rat_from_human():
    assert uniprot_client.get_rat_id('P04049') == 'P11345'


def test_length():
    assert uniprot_client.get_length('P15056') == 766


@attr('webservice')
def test_get_function():
    fun = uniprot_client.get_function('P15056')
    assert fun.startswith('Protein kinase involved in the transduction')


def test_get_is_reviewed():
    assert not uniprot_client.is_reviewed('V9HWD6')
    assert uniprot_client.is_reviewed('P31946-1')


def test_get_ids_from_refseq():
    up_ids = uniprot_client.get_ids_from_refseq('NP_003395.1')
    assert set(up_ids) == {'V9HWD6', 'P31946-1'}, set(up_ids)
    up_ids = uniprot_client.get_ids_from_refseq('NP_003395.1',
                                                reviewed_only=True)
    assert up_ids == ['P31946-1']


@attr('webservice')
def test_get_signal_peptide():
    # This is a valid entry local to the resource file
    sp = uniprot_client.get_signal_peptide('P00533')
    assert sp.begin == 1, sp
    assert sp.end == 24, sp
    # This one requires a web lookup
    sp = uniprot_client.get_signal_peptide('P00534')
    assert sp is None, sp
    # This one errors when doing web lookup
    sp = uniprot_client.get_signal_peptide('Q9H7H1')
    assert sp is None, sp


def test_get_hgnc_id():
    hgnc_id = uniprot_client.get_hgnc_id('P07305')
    assert hgnc_id == '4714', hgnc_id
    # NRXN2: ['P58401', 'Q9P2S2'] 8009
    hgnc_id = uniprot_client.get_hgnc_id('P58401')
    assert hgnc_id == '8009', hgnc_id
    hgnc_id = uniprot_client.get_hgnc_id('Q9P2S2')
    assert hgnc_id == '8009', hgnc_id


def test_is_mouse():
    assert uniprot_client.is_mouse('P28028-1') is True
    assert uniprot_client.is_mouse('P07305') is False


def test_is_rat():
    assert uniprot_client.is_rat('P11345-1') is True
    assert uniprot_client.is_rat('P28028') is False


def test_process_chain():
    chain_str = ('CHAIN 1..7096;  /note="Replicase polyprotein 1ab";  '
                 '/id="PRO_0000449618";  CHAIN 1..180;  /note="Host '
                 'translation inhibitor nsp1";  /id="PRO_0000449619";')
    chains = _process_feature('CHAIN', chain_str)
    assert len(chains) == 2
    assert chains[0].id == 'PRO_0000449618'
    assert chains[0].begin == 1, chains
    assert chains[0].end == 7096
    assert chains[0].name == 'Replicase polyprotein 1ab'
    assert chains[1].id == 'PRO_0000449619', chains
    assert chains[1].begin == 1
    assert chains[1].end == 180
    assert chains[1].name == 'Host translation inhibitor nsp1'


def test_features():
    features = uniprot_client.get_features('P55957')
    assert len(features) == 4, features
    chains = uniprot_client.get_chains('P55957')
    assert len(chains) == 4
    assert 'BH3-interacting domain death agonist p15' in \
           {c.name for c in chains}
    for chain in chains:
        assert chain.type == 'CHAIN'
        if chain.name == 'BH3-interacting domain death agonist p15':
            assert chain.begin == 62, chain
            assert chain.end == 195
            assert chain.id == 'PRO_0000223233'


def test_feature_by_id():
    feature = uniprot_client.get_feature_by_id('PRO_0000292268')
    assert feature is not None
    assert feature.name == 'Processed angiotensin-converting enzyme 2'
    assert feature.type == 'CHAIN'
    assert feature.begin == 18, feature
    assert feature.end == 708


def test_sars_cov2_feature():
    feat = uniprot_client.get_feature_by_id('PRO_0000449635')
    assert feat.type == 'CHAIN'
    assert feat.begin == 1, feat
    assert feat.end == 180, feat
    assert feat.name == 'Host translation inhibitor nsp1', feat.name


def test_get_feature_of():
    up_id = uniprot_client.get_feature_of('PRO_0000449635')
    assert up_id == 'P0DTC1', up_id


def test_entrez_uniprot():
    assert uniprot_client.get_entrez_id('Q66K41') == '201181'
    assert uniprot_client.get_id_from_entrez('201181') == 'Q66K41'
