import re
import csv
import gzip
import json
import logging
import itertools
import requests
from typing import List, Tuple, Optional, Union
from functools import lru_cache
import xml
from xml.etree import ElementTree
from urllib.error import HTTPError
from protmapper.resources import resource_manager, feature_from_json, Feature

logger = logging.getLogger(__name__)


uniprot_url = 'https://legacy.uniprot.org/uniprot/'

xml_ns = {'up': 'http://uniprot.org/uniprot'}


@lru_cache(maxsize=10000)
def query_protein(protein_id: str) -> Union[ElementTree.ElementTree, None]:
    """Retrieve the XML entry for a given protein.

    Parameters
    ----------
    protein_id :
        The UniProt ID of the protein to look up.

    Returns
    -------
    :
        An ElementTree representation of the XML entry for the
        protein.
    """
    # Try looking up a primary ID if the given one
    # is a secondary ID and strip off isoforms
    protein_id = get_primary_id(_strip_isoform(protein_id))
    url = uniprot_url + protein_id + '.xml'
    try:
        # As opposed to the RDF endpoint, the XML endpoint returns
        # an identical entry for secondary accessions, for instance,
        # the response for the secondary ID A0A021WW06 is identical to
        # the response for the primary ID P40417.
        ret = requests.get(url)
        et = ElementTree.fromstring(ret.content)
        return et
    except Exception as e:
        return None


def _strip_isoform(protein_id):
    return protein_id.split('-')[0]


def _split_isoform(protein_id):
    parts = protein_id.split('-', maxsplit=1)
    protein_id = parts[0]
    isoform = None
    if len(parts) == 2:
        if re.match(r'\d+', parts[1]):
            isoform = parts[1]
    return protein_id, isoform


def _reattach_isoform(pid, iso):
    if iso is not None:
        return '%s-%s' % (pid, iso)
    else:
        return pid


def is_secondary(protein_id):
    """Return True if the UniProt ID corresponds to a secondary accession.

    Parameters
    ----------
    protein_id : str
        The UniProt ID to check.

    Returns
    -------
    True if it is a secondary accessing entry, False otherwise.
    """
    entry = um.uniprot_sec.get(_strip_isoform(protein_id))
    if not entry:
        return False
    return True


def is_reviewed(protein_id):
    """Return True if the UniProt ID corresponds to a reviewed entry.

    Parameters
    ----------
    protein_id : str
        The UniProt ID to check.

    Returns
    -------
    True if it is a reviewed entry, False otherwise.
    """
    return _strip_isoform(protein_id) in um.uniprot_reviewed


def get_primary_id(protein_id):
    """Return a primary entry corresponding to the UniProt ID.

    Parameters
    ----------
    protein_id : str
        The UniProt ID to map to primary.

    Returns
    -------
    primary_id : str
        If the given ID is primary, it is returned as is. Otherwise the primary
        IDs are looked up. If there are multiple primary IDs then the first
        human one is returned. If there are no human primary IDs then the
        first primary found is returned.
    """
    base_id, isoform = _split_isoform(protein_id)
    primaries = um.uniprot_sec.get(base_id)
    if primaries:
        if len(primaries) > 1:
            logger.debug('More than 1 primary ID for %s.' % base_id)
            for primary in primaries:
                # Often secondary IDs were broken into multiple primary IDs
                # for different organisms. In this case we return the human
                # one if it exists.
                if is_human(primary):
                    return _reattach_isoform(primary, isoform)
        # If we haven't returned anything then we just return the
        # first primary id
        return _reattach_isoform(primaries[0], isoform)
    # If there is no secondary entry then we assume this is a primary entry
    return protein_id


def get_family_members(family_name, human_only=True):
    """Return the HGNC gene symbols which are the members of a given family.

    Parameters
    ----------
    family_name : str
        Family name to be queried.
    human_only : bool
        If True, only human proteins in the family will be returned.
        Default: True

    Returns
    -------
    gene_names : list
        The HGNC gene symbols corresponding to the given family.
    """
    data = {'query': 'family:%s' % family_name,
            'format': 'list'}
    if human_only:
        data['fil'] = 'organism:human'
    res = requests.get(uniprot_url, params=data)
    if not res.status_code == 200 or not res.text:
        return None
    # res.text gets us the Unicode
    html = res.text
    protein_list = html.strip().split('\n')
    gene_names = []
    for p in protein_list:
        gene_name = get_gene_name(p)
        gene_names.append(gene_name)
    return gene_names


def get_mnemonic(protein_id, web_fallback=False):
    """Return the UniProt mnemonic for the given UniProt ID.

    Parameters
    ----------
    protein_id : str
        UniProt ID to be mapped.
    web_fallback : Optional[bool]
        If True and the offline lookup fails, the UniProt web service
        is used to do the query.

    Returns
    -------
    mnemonic : str
        The UniProt mnemonic corresponding to the given Uniprot ID.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    try:
        mnemonic = um.uniprot_mnemonic[protein_id]
        return mnemonic
    except KeyError:
        pass
    if not web_fallback:
        return None

    tree = query_protein(protein_id)
    if tree is None:
        return None

    mnemonic = tree.find('up:entry/up:name', namespaces=xml_ns)
    if mnemonic is None:
        return None
    return mnemonic.text


def get_id_from_mnemonic(uniprot_mnemonic):
    """Return the UniProt ID for the given UniProt mnemonic.

    Parameters
    ----------
    uniprot_mnemonic : str
        UniProt mnemonic to be mapped.

    Returns
    -------
    uniprot_id : str
        The UniProt ID corresponding to the given Uniprot mnemonic.
    """
    try:
        uniprot_id = um.uniprot_mnemonic_reverse[uniprot_mnemonic]
        return uniprot_id
    except KeyError:
        return None


def get_gene_name(protein_id, web_fallback=True):
    """Return the gene name or canonical protein name for the given UniProt ID.

    If available, this function returns the primary gene name provided by
    UniProt. If not available, the primary protein name is returned.

    Parameters
    ----------
    protein_id : str
        UniProt ID to be mapped.
    web_fallback : Optional[bool]
        If True and the offline lookup fails, the UniProt web service
        is used to do the query.

    Returns
    -------
    gene_name : str
        The gene name corresponding to the given Uniprot ID.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    try:
        gene_name = um.uniprot_gene_name[protein_id]
        # We only get here if the protein_id was in the dict
        if gene_name:
            return gene_name
        # We do it this way to return None for empty strings
        else:
            return None
    except KeyError:
        if not web_fallback:
            return None

    tree = query_protein(protein_id)
    if tree is None:
        return None
    name = tree.find('up:entry/up:gene/up:name', namespaces=xml_ns)
    if name is not None:
        return name.text
    return None


def get_gene_synonyms(protein_id: str) -> List[str]:
    """Return a list of synonyms for the gene corresponding to a protein.

    Note that synonyms here also include the official gene name as
    returned by get_gene_name.

    Parameters
    ----------
    protein_id :
        The UniProt ID of the protein to query

    Returns
    -------
    :
        The list of synonyms of the gene corresponding to the protein
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    protein = query_protein(protein_id)
    if protein is None:
        return []
    synonyms = []
    gene_synoyms = protein.findall('up:entry/up:gene/up:name',
                                   namespaces=xml_ns)
    for gene_syn in gene_synoyms:
        synonyms.append(gene_syn.text)
    return synonyms


def get_protein_synonyms(protein_id):
    """Return a list of synonyms for a protein.

    Note that this function returns protein synonyms as provided by UniProt.
    The get_gene_synonym returns synonyms given for the gene corresponding
    to the protein, and get_synonyms returns both.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein to query

    Returns
    -------
    synonyms : list[str]
        The list of synonyms of the protein
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    tree = query_protein(protein_id)
    if tree is None:
        return None
    synonyms = []
    for syn_type, syn_len in itertools.product(['recommended', 'alternative'],
                                               ['full', 'short']):
        synonym_type = 'up:entry/up:protein/up:%sName/up:%sName' % \
            (syn_type, syn_len)
        synonyms_xml = tree.findall(synonym_type, namespaces=xml_ns)
        for synonym_xml in synonyms_xml:
            synonyms.append(synonym_xml.text)
    return synonyms


def get_synonyms(protein_id):
    """Return synonyms for a protein and its associated gene.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein to query

    Returns
    -------
    synonyms : list[str]
        The list of synonyms of the protein and its associated gene.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    ret = []
    gene_syms = get_gene_synonyms(protein_id)
    if gene_syms:
        ret.extend(gene_syms)
    prot_syms = get_protein_synonyms(protein_id)
    if prot_syms:
        ret.extend(prot_syms)
    return ret


@lru_cache(maxsize=10000)
def get_sequence(protein_id):
    base, iso = _split_isoform(get_primary_id(protein_id))
    # Try to get the sequence from the downloaded sequence files
    if iso == '1':
        protein_id = base
    else:
        protein_id = _reattach_isoform(base, iso)
    seq = um.uniprot_sequences.get(protein_id)
    if seq is None:
        url = uniprot_url + '%s.fasta' % protein_id
        res = requests.get(url)
        res.raise_for_status()
        # res.text is Unicode
        lines = res.text.splitlines()
        seq = (''.join(lines[1:])).replace('\n', '')
    return seq


def get_modifications(protein_id: str) -> List[Tuple[str, int]]:
    """Return a list of modifications for a protein.

    Parameters
    ----------
    protein_id :
        The UniProt ID of the protein to query

    Returns
    -------
    :
        The list of modifications of the protein, each represented
        as a tuple of residue description string and position
        string.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    tree = query_protein(protein_id)
    if tree is None:
        return None

    # We find all features of type 'modified residue'
    features = tree.findall("up:entry/up:feature[@type='modified residue']",
                            namespaces=xml_ns)
    mods = []
    for feature in features:
        # We find the position of the modified residue
        pos_tag = feature.find('up:location/up:position', namespaces=xml_ns)
        if pos_tag is None:
            continue
        pos = int(pos_tag.attrib['position'])
        # We find the residue
        res = feature.attrib['description'].split(';')[0]
        mods.append((res, pos))

    return mods


def verify_location(protein_id, residue, location):
    """Return True if the residue is at the given location in the UP sequence.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein whose sequence is used as reference.
    residue : str
        A single character amino acid symbol (Y, S, T, V, etc.)
    location : str
        The location on the protein sequence (starting at 1) at which the
        residue should be checked against the reference sequence.

    Returns
    -------
    True if the given residue is at the given position in the sequence
    corresponding to the given UniProt ID, otherwise False.
    """
    seq = get_sequence(protein_id)
    # If we couldn't get the sequence (can happen due to web service hiccups)
    # don't throw the statement away by default
    if seq is None:
        return True
    try:
        loc_int = int(location)
    except ValueError:
        logger.warning('Invalid location %s' % location)
        loc_int = -1

    if (loc_int < 1) or (loc_int > len(seq)):
        return False
    elif seq[loc_int - 1] == residue:
        return True
    return False


def verify_modification(protein_id, residue, location=None):
    """Return True if the residue at the given location has a known modifiation.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein whose sequence is used as reference.
    residue : str
        A single character amino acid symbol (Y, S, T, V, etc.)
    location : Optional[str]
        The location on the protein sequence (starting at 1) at which the
        modification is checked.

    Returns
    -------
    True if the given residue is reported to be modified at the given position
    in the sequence corresponding to the given UniProt ID, otherwise False.
    If location is not given, we only check if there is any residue of the
    given type that is modified.
    """
    mods = get_modifications(protein_id)
    mod_locs = [m[1] for m in mods]
    if location:
        if not verify_location(protein_id, residue, location):
            return False
        try:
            mod_idx = mod_locs.index(location)
        except ValueError:
            return False
        return True
    else:
        seq = get_sequence(protein_id)
        for ml in mod_locs:
            if seq[ml - 1] == residue:
                return True
        return False


def _is_organism(protein_id, organism_suffix):
    protein_id = get_primary_id(_strip_isoform(protein_id))
    mnemonic = get_mnemonic(protein_id)
    if mnemonic is None:
        return False
    if mnemonic.endswith(organism_suffix):
        return True
    return False


def is_human(protein_id):
    """Return True if the given protein id corresponds to a human protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein

    Returns
    -------
    True if the protein_id corresponds to a human protein, otherwise False.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return _is_organism(protein_id, 'HUMAN')


def is_mouse(protein_id):
    """Return True if the given protein id corresponds to a mouse protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein

    Returns
    -------
    True if the protein_id corresponds to a mouse protein, otherwise False.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return _is_organism(protein_id, 'MOUSE')


def is_rat(protein_id):
    """Return True if the given protein id corresponds to a rat protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein

    Returns
    -------
    True if the protein_id corresponds to a rat protein, otherwise False.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return _is_organism(protein_id, 'RAT')


def get_hgnc_id(protein_id):
    """Return the HGNC ID given the protein id of a human protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the human protein

    Returns
    -------
    hgnc_id : str
        HGNC ID of the human protein
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return um.uniprot_hgnc.get(protein_id)


def get_entrez_id(protein_id):
    """Return the Entrez ID given a protein ID.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the protein

    Returns
    -------
    str or None
        Entrez ID of the corresponding gene or None if not available.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return um.uniprot_entrez.get(protein_id)


def get_id_from_entrez(entrez_id):
    """Return the UniProt ID given the Entrez ID of a gene.

    Parameters
    ----------
    entrez_id : str
        Entrez ID of the gene

    Returns
    -------
    str or None
        UniProt ID of the corresponding protein or None if not available.
    """
    return um.entrez_uniprot.get(entrez_id)


def get_mgi_id(protein_id):
    """Return the MGI ID given the protein id of a mouse protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the mouse protein

    Returns
    -------
    mgi_id : str
        MGI ID of the mouse protein
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return um.uniprot_mgi.get(protein_id)


def get_rgd_id(protein_id):
    """Return the RGD ID given the protein id of a rat protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of the rat protein

    Returns
    -------
    rgd_id : str
        RGD ID of the rat protein
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return um.uniprot_rgd.get(protein_id)


def get_id_from_mgi(mgi_id):
    """Return the UniProt ID given the MGI ID of a mouse protein.

    Parameters
    ----------
    mgi_id : str
        The MGI ID of the mouse protein.

    Returns
    -------
    up_id : str
        The UniProt ID of the mouse protein.
    """
    return um.uniprot_mgi_reverse.get(mgi_id)


def get_id_from_mgi_name(mgi_name: str) -> Optional[str]:
    """Return the UniProt ID given the MGI name of a mouse protein.

    Parameters
    ----------
    mgi_name : str
        The MGI name of the mouse protein.

    Returns
    -------
    up_id : str
        The UniProt ID of the mouse protein.
    """
    return um.mgi_name_to_up.get(mgi_name)


def get_id_from_rgd(rgd_id):
    """Return the UniProt ID given the RGD ID of a rat protein.

    Parameters
    ----------
    rgd_id : str
        The RGD ID of the rat protein.

    Returns
    -------
    up_id : str
        The UniProt ID of the rat protein.
    """
    return um.uniprot_rgd_reverse.get(rgd_id)


def get_id_from_rgd_name(rgd_name: str) -> Optional[str]:
    """Return the UniProt ID given the RGD name of a rat protein.

    Parameters
    ----------
    rgd_name : str
        The RGD name of the rat protein.

    Returns
    -------
    up_id : str
        The UniProt ID of the rat protein.
    """
    return um.rgd_name_to_up.get(rgd_name)


def get_mouse_id(human_protein_id):
    """Return the mouse UniProt ID given a human UniProt ID.

    Parameters
    ----------
    human_protein_id : str
        The UniProt ID of a human protein.

    Returns
    -------
    mouse_protein_id : str
        The UniProt ID of a mouse protein orthologous to the given human
        protein.
    """
    human_protein_id = get_primary_id(_strip_isoform(human_protein_id))
    return um.uniprot_human_mouse.get(human_protein_id)


def get_rat_id(human_protein_id):
    """Return the rat UniProt ID given a human UniProt ID.

    Parameters
    ----------
    human_protein_id : str
        The UniProt ID of a human protein.

    Returns
    -------
    rat_protein_id : str
        The UniProt ID of a rat protein orthologous to the given human protein
    """
    human_protein_id = get_primary_id(_strip_isoform(human_protein_id))
    return um.uniprot_human_rat.get(human_protein_id)


def get_length(protein_id):
    """Return the length (number of amino acids) of a protein.

    Parameters
    ----------
    protein_id : str
        UniProt ID of a protein.

    Returns
    -------
    length : int
        The length of the protein in amino acids.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    return um.uniprot_length.get(protein_id)


def get_function(protein_id):
    """Return the function description of a given protein.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein.

    Returns
    -------
    str
        The function description of the protein.
    """
    et = query_protein(protein_id)
    if et is None:
        return None
    function = et.find('up:entry/up:comment[@type="function"]/up:text',
                       namespaces=xml_ns)
    if function is None:
        return None
    return function.text


def get_features(protein_id):
    """Return a list of features (chains, peptides) for a given protein.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein whose features are to be returned.

    Returns
    -------
    list of Feature
        A list of Feature named tuples representing each Feature.
    """
    return um.features.get(protein_id, [])


def get_chains(protein_id):
    """Return the list of cleaved chains for the given protein.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein whose cleaved chains are to be returned.

    Returns
    -------
    list of Feature
        A list of Feature named tuples representing each chain.
    """
    features = get_features(protein_id)
    chains = [f for f in features if f.type == 'CHAIN']
    return chains


def get_signal_peptide(protein_id, web_fallback=True):
    """Return the position of a signal peptide for the given protein.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein whose signal peptide position
        is to be returned.
    web_fallback : Optional[bool]
        If True the UniProt web service is used to download information when
        the local resource file doesn't contain the right information.

    Returns
    -------
    Feature
        A Feature named tuple representing the signal peptide.
    """
    protein_id = get_primary_id(_strip_isoform(protein_id))
    # Note, we use False here to differentiate from None
    if not web_fallback and protein_id not in um.features:
        return None
    elif protein_id in um.features:
        sp = [f for f in get_features(protein_id) if f.type == 'SIGNAL']
        if sp:
            return sp[0]
        elif not web_fallback:
            return False

    et = query_protein(protein_id)
    if et is None:
        return None
    location = et.find(
        'up:entry/up:feature[@type="signal peptide"]/up:location',
        namespaces=xml_ns)
    begin_pos = None
    end_pos = None
    if location is not None:
        begin = location.find('up:begin', namespaces=xml_ns)
        if begin is not None:
            begin_pos = begin.attrib.get('position')
            if begin_pos is not None:
                begin_pos = int(begin_pos)
        end = location.find('up:end', namespaces=xml_ns)
        if end is not None:
            end_pos = end.attrib.get('position')
            if end_pos is not None:
                end_pos = int(end_pos)
    if begin_pos is not None and end_pos is not None:
        return Feature('SIGNAL', begin_pos, end_pos, None, None)
    return None


def get_feature_by_id(feature_id):
    """Return a Feature based on its unique feature ID.

    Parameters
    ----------
    feature_id : str
        A Feature ID, of the form PRO_*.

    Returns
    -------
    Feature or None
        A Feature with the given ID.
    """
    return um.features_by_id.get(feature_id)


def get_feature_of(feature_id):
    """Return the UniProt ID of the protein to which the given feature belongs.

    Parameters
    ----------
    feature_id : str
        A Feature ID, of the form PRO_*.

    Returns
    -------
    str or None
        A UniProt ID corresponding to the given feature, or None if not
        available (generally shouldn't happen, unless the feature ID is
        invalid).
    """
    for up_id, feats in um.features.items():
        for feat in feats:
            if feat.id == feature_id:
                return up_id
    return None


def get_ids_from_refseq(refseq_id, reviewed_only=False):
    """Return UniProt IDs from a RefSeq ID".

    Parameters
    ----------
    refseq_id : str
        The RefSeq ID of the protein to map.
    reviewed_only : Optional[bool]
        If True, only reviewed UniProt IDs are returned.
        Default: False

    Returns
    -------
    list of str
        A list of UniProt IDs corresponding to the RefSeq ID.
    """
    try:
        up_ids = um.refseq_uniprot[refseq_id]
    except KeyError:
        return []
    primaries = list(set([get_primary_id(up_id) for up_id in up_ids]))
    if reviewed_only:
        return [up_id for up_id in primaries if is_reviewed(up_id)]
    else:
        return primaries


def get_organism_id(protein_id):
    """Return the Taxonomy ID of the organism that a protein belongs to.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of a protein.

    Returns
    -------
    str or None
        The Taxonomy ID of the organism the protein belongs to or None
        if not available.
    """
    return um.organism_ids.get(protein_id)


class UniprotMapper(object):
    def __init__(self):
        self.initialized = False
        self.initialized_seq = False
        self.initialized_hgnc = False
        self.initialized_refseq = False
        self._entrez_uniprot = {}
        self._uniprot_entrez = {}

    def initialize(self):
        maps = _build_uniprot_entries()
        (self._uniprot_gene_name, self._uniprot_mnemonic,
         self._uniprot_mnemonic_reverse, self._uniprot_mgi,
         self._uniprot_rgd, self._uniprot_mgi_reverse,
         self._uniprot_rgd_reverse, self._uniprot_length,
         self._uniprot_reviewed, self._features, self._features_by_id,
         self._organisms_by_id, _uniprot_entrez,
         _entrez_uniprot, self._mgi_name_to_up, self._rgd_name_to_up) = maps

        # Here we don't overwrite the value to de-prioritize
        for k, v in _uniprot_entrez.items():
            if k not in self._uniprot_entrez:
                self._uniprot_entrez[k] = v
        for k, v in _entrez_uniprot.items():
            if k not in self._entrez_uniprot:
                self._entrez_uniprot[k] = v

        self._uniprot_sec = _build_uniprot_sec()

        self.initialized = True

    def initialize_hgnc(self):
        self._uniprot_human_mouse, self._uniprot_human_rat = \
            _build_human_mouse_rat()
        _, _, self._uniprot_hgnc, _entrez_uniprot, \
            _uniprot_entrez = _build_hgnc_mappings()

        # Here we always overwrite the value to prioritize
        for upid, egid in _uniprot_entrez.items():
            for uid in upid.split(', '):
                self._uniprot_entrez[uid] = egid
        for egid, upid in _entrez_uniprot.items():
            if ',' not in upid:
                self._entrez_uniprot[egid] = upid

        self.initialized_hgnc = True

    def initialize_seq(self):
        self._sequences = _build_uniprot_sequences()
        self.initialized_seq = True

    def initialize_refseq(self):
        self._refseq_uniprot = _build_refseq_uniprot()
        self.initialized_refseq = True

    @property
    def uniprot_gene_name(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_gene_name

    @property
    def uniprot_hgnc(self):
        if not self.initialized_hgnc:
            self.initialize_hgnc()
        return self._uniprot_hgnc

    @property
    def uniprot_entrez(self):
        if not self.initialized:
            self.initialize()
        if not self.initialized_hgnc:
            self.initialize_hgnc()
        return self._uniprot_entrez

    @property
    def entrez_uniprot(self):
        if not self.initialized:
            self.initialize()
        if not self.initialized_hgnc:
            self.initialize_hgnc()
        return self._entrez_uniprot

    @property
    def uniprot_mnemonic(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_mnemonic

    @property
    def uniprot_mnemonic_reverse(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_mnemonic_reverse

    @property
    def uniprot_mgi(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_mgi

    @property
    def uniprot_rgd(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_rgd

    @property
    def uniprot_mgi_reverse(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_mgi_reverse

    @property
    def uniprot_rgd_reverse(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_rgd_reverse

    @property
    def uniprot_length(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_length

    @property
    def uniprot_reviewed(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_reviewed

    @property
    def uniprot_sec(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_sec

    @property
    def uniprot_human_mouse(self):
        if not self.initialized_hgnc:
            self.initialize_hgnc()
        return self._uniprot_human_mouse

    @property
    def uniprot_human_rat(self):
        if not self.initialized_hgnc:
            self.initialize_hgnc()
        return self._uniprot_human_rat

    @property
    def uniprot_sequences(self):
        if not self.initialized_seq:
            self.initialize_seq()
        return self._sequences

    @property
    def refseq_uniprot(self):
        if not self.initialized_refseq:
            self.initialize_refseq()
        return self._refseq_uniprot

    @property
    def features(self):
        if not self.initialized:
            self.initialize()
        return self._features

    @property
    def features_by_id(self):
        if not self.initialized:
            self.initialize()
        return self._features_by_id

    @property
    def organism_ids(self):
        if not self.initialized:
            self.initialize()
        return self._organisms_by_id

    @property
    def mgi_name_to_up(self):
        if not self.initialized:
            self.initialize()
        return self._mgi_name_to_up

    @property
    def rgd_name_to_up(self):
        if not self.initialized:
            self.initialize()
        return self._rgd_name_to_up


um = UniprotMapper()


def _build_uniprot_entries():
    up_entries_file = resource_manager.get_create_resource_file('up')
    uniprot_gene_name = {}
    uniprot_mnemonic = {}
    uniprot_mnemonic_reverse = {}
    uniprot_mgi = {}
    uniprot_rgd = {}
    uniprot_mgi_reverse = {}
    uniprot_rgd_reverse = {}
    uniprot_length = {}
    uniprot_features = {}
    uniprot_reviewed = set()
    organisms_by_id = {}
    uniprot_entrez = {}
    uniprot_entrez_reverse = {}
    files = [up_entries_file]
    mgi_name_to_up = {}
    rgd_name_to_up = {}
    for file in files:
        with gzip.open(file, 'rt', encoding='utf-8') as fh:
            csv_rows = csv.reader(fh, delimiter='\t')
            # Skip the header row
            next(csv_rows)
            for row in csv_rows:
                up_id, gene_name, up_mnemonic, rgd, mgi, length, reviewed, \
                    organism_id, entrez_id, features_json = row
                # Store the entry in the reviewed set
                if reviewed == 'reviewed':
                    uniprot_reviewed.add(up_id)
                # This is to turn empty strings into explicit Nones
                uniprot_gene_name[up_id] = gene_name if gene_name else None
                uniprot_mnemonic[up_id] = up_mnemonic
                uniprot_mnemonic_reverse[up_mnemonic] = up_id
                uniprot_length[up_id] = int(length)
                if mgi:
                    mgi_ids = mgi.split(';')
                    if mgi_ids:
                        uniprot_mgi[up_id] = mgi_ids[0]
                        uniprot_mgi_reverse[mgi_ids[0]] = up_id
                        mgi_name_to_up[gene_name] = up_id
                if rgd:
                    rgd_ids = rgd.split(';')
                    if rgd_ids:
                        uniprot_rgd[up_id] = rgd_ids[0]
                        uniprot_rgd_reverse[rgd_ids[0]] = up_id
                        rgd_name_to_up[gene_name] = up_id
                uniprot_features[up_id] = [feature_from_json(feat) for
                                           feat in json.loads(features_json)]
                organisms_by_id[up_id] = organism_id

                # Entrez mappings
                entrez_ids = [ei for ei in
                              [e.strip() for e in entrez_id.split(';')] if ei]
                for eid in entrez_ids:
                    uniprot_entrez[up_id] = eid
                    uniprot_entrez_reverse[eid] = up_id

    # Build a dict of features by feature ID
    features_by_id = {}
    for up_id, feats in uniprot_features.items():
        for feat in feats:
            features_by_id[feat.id] = feat

    return (uniprot_gene_name, uniprot_mnemonic, uniprot_mnemonic_reverse,
            uniprot_mgi, uniprot_rgd, uniprot_mgi_reverse, uniprot_rgd_reverse,
            uniprot_length, uniprot_reviewed, uniprot_features, features_by_id,
            organisms_by_id, uniprot_entrez, uniprot_entrez_reverse,
            mgi_name_to_up, rgd_name_to_up)


def _build_human_mouse_rat():
    hgnc_file = resource_manager.get_create_resource_file('hgnc')
    with gzip.open(hgnc_file, 'rt', encoding='utf-8') as fh:
        csv_rows = csv.reader(fh, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        uniprot_mouse = {}
        uniprot_rat = {}
        for row in csv_rows:
            human_id, mgi_id, rgd_id = row[6:9]
            if human_id:
                if mgi_id:
                    mgi_id = mgi_id.split(', ')[0]
                    if mgi_id.startswith('MGI:'):
                        mgi_id = mgi_id[4:]
                    mouse_id = um.uniprot_mgi_reverse.get(mgi_id)
                    if mouse_id:
                        uniprot_mouse[human_id] = mouse_id
                if rgd_id:
                    rgd_id = rgd_id.split(', ')[0]
                    if rgd_id.startswith('RGD:'):
                        rgd_id = rgd_id[4:]
                    rat_id = um.uniprot_rgd_reverse.get(rgd_id)
                    if rat_id:
                        uniprot_rat[human_id] = rat_id
    return uniprot_mouse, uniprot_rat


def _build_hgnc_mappings():
    hgnc_file = resource_manager.get_create_resource_file('hgnc')
    with gzip.open(hgnc_file, 'rt', encoding='utf-8') as fh:
        csv_rows = csv.reader(fh, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        hgnc_name_to_id = {}
        hgnc_id_to_up = {}
        up_to_hgnc_id = {}
        entrez_to_up = {}
        up_to_entrez = {}
        for row in csv_rows:
            hgnc_id = row[0][5:]
            hgnc_status = row[3]
            if hgnc_status == 'Approved':
                hgnc_name = row[1]
                hgnc_name_to_id[hgnc_name] = hgnc_id
            # Uniprot
            uniprot_id = row[6]
            if uniprot_id:
                hgnc_id_to_up[hgnc_id] = uniprot_id
                uniprot_ids = uniprot_id.split(', ')
                for upid in uniprot_ids:
                    up_to_hgnc_id[upid] = hgnc_id
                # Entrez
                entrez_id = row[5]
                if entrez_id:
                    for upid in uniprot_ids:
                        up_to_entrez[upid] = entrez_id
                        entrez_to_up[entrez_id] = uniprot_id

    return hgnc_name_to_id, hgnc_id_to_up, up_to_hgnc_id, \
        entrez_to_up, up_to_entrez


def _build_uniprot_sec():
    # File containing secondary accession numbers mapped
    # to primary accession numbers
    sec_file = resource_manager.get_create_resource_file('upsec')
    uniprot_sec = {}
    with gzip.open(sec_file, 'rt', encoding='utf-8') as fh:
        lines = fh.readlines()
        for i, l in enumerate(lines):
            if l.startswith('Secondary AC'):
                entry_lines = lines[i+2:]

    for l in entry_lines:
        sec_id, prim_id = l.split()
        try:
            uniprot_sec[sec_id].append(prim_id)
        except KeyError:
            uniprot_sec[sec_id] = [prim_id]
    return uniprot_sec


def _build_uniprot_sequences():
    seq_file = resource_manager.get_create_resource_file('swissprot',
                                                         cached=True)
    iso_file = resource_manager.get_create_resource_file('isoforms',
                                                         cached=True)
    logger.info("Loading Swissprot sequences...")
    sp_seq = load_fasta_sequences(seq_file)
    logger.info("Loading Uniprot isoform sequences...")
    iso_seq = load_fasta_sequences(iso_file)
    sp_seq.update(iso_seq)
    return sp_seq


def _build_refseq_uniprot():
    refseq_uniprot_file = resource_manager.get_create_resource_file(
                                                'refseq_uniprot')
    refseq_up = {}
    with gzip.open(refseq_uniprot_file, 'rt', encoding='utf-8') as f:
        csvreader = csv.reader(f)
        for refseq_id, up_id in csvreader:
            if refseq_id not in refseq_up:
                refseq_up[refseq_id] = []
            refseq_up[refseq_id].append(up_id)
    return refseq_up


def load_fasta_sequences(seq_file, id_delimiter='|', id_index=1):
    if seq_file.endswith('gz'):
        with gzip.open(seq_file, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
    # This is necessary for downstream usage of this function
    else:
        with open(seq_file, 'rt', encoding='utf-8') as f:
            lines = f.readlines()
    return load_fasta_sequence_lines(lines, id_delimiter=id_delimiter,
                                     id_index=id_index)


def load_fasta_sequence_lines(lines, id_delimiter='|', id_index=1):
    sequences = {}
    cur_id = None
    seq_lines = []
    for line in lines:
        if line.startswith('>'):
            line_id = line[1:].split(id_delimiter)[id_index]
            if cur_id is not None:
                seq = ''.join(seq_lines)
                sequences[cur_id] = seq
                seq_lines = []
            cur_id = line_id
        else:
            seq_lines.append(line.strip())
    # Add the last sequence
    seq = ''.join(seq_lines)
    sequences[cur_id] = seq
    return sequences
