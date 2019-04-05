import os
import re
import csv
import rdflib
import logging
import requests
from xml.etree import ElementTree
from functools import lru_cache
from urllib.error import HTTPError
from protmapper.resources import resource_manager

logger = logging.getLogger(__name__)


uniprot_url = 'http://www.uniprot.org/uniprot/'


rdf_prefixes = """
    PREFIX up: <http://purl.uniprot.org/core/>
    PREFIX db: <http://purl.uniprot.org/database/>
    PREFIX faldo: <http://biohackathon.org/resource/faldo#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#> """

xml_ns = {'up': 'http://uniprot.org/uniprot'}


@lru_cache(maxsize=10000)
def query_protein(protein_id):
    """Return the UniProt entry as an RDF graph for the given UniProt ID.

    Parameters
    ----------
    protein_id : str
        UniProt ID to be queried.

    Returns
    -------
    g : rdflib.Graph
        The RDF graph corresponding to the UniProt entry.
    """
    # Try looking up a primary ID if the given one
    # is a secondary ID
    try:
        prim_ids = um.uniprot_sec[protein_id]
        protein_id = prim_ids[0]
    except KeyError:
        pass
    url = uniprot_url + protein_id + '.rdf'
    g = rdflib.Graph()
    try:
        g.parse(url)
    except HTTPError:
        logger.warning('Could not find protein with id %s' % protein_id)
        return None
    except rdflib.exceptions.ParserError as e:
        logger.error('Could not parse RDF at %s' % url)
        logger.error(e)
        return None

    # Check if the entry has been replaced by a new entry
    query = rdf_prefixes + """
        SELECT ?res2
        WHERE {
            ?res1 up:replacedBy ?res2 .
            }
        """
    res = g.query(query)
    if res:
        term = [r for r in res][0][0]
        replaced_by_id = term.split('/')[-1]
        return query_protein(replaced_by_id)
    return g


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
    entry = um.uniprot_sec.get(protein_id)
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
    # If this is an isoform, check to see if the base ID is reviewed
    # (isoform IDs are not in the ID list)
    if '-' in protein_id:
        protein_id = protein_id.split('-')[0]
    return (protein_id in um.uniprot_reviewed)


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
    primaries = um.uniprot_sec.get(protein_id)
    if primaries:
        if len(primaries) > 1:
            logger.debug('More than 1 primary ID for %s.' % protein_id)
            for primary in primaries:
                # Often secondary IDs were broken into multiple primary IDs
                # for different organisms. In this case we return the human
                # one if it exists.
                if is_human(primary):
                    return primary
        # If we haven't returned anything then we just return the
        # first primary id
        return primaries[0]
    # If there is not secondary entry the we assume this is a primary entry
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
    try:
        mnemonic = um.uniprot_mnemonic[protein_id]
        return mnemonic
    except KeyError:
        pass
    if not web_fallback:
        return None
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?mnemonic
        WHERE {
            ?r up:mnemonic ?mnemonic .
        }
        """
    res = g.query(query)
    if res:
        mnemonic = [r for r in res][0][0].toPython()
        return mnemonic
    else:
        return None


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
    """Return the gene name for the given UniProt ID.

    This is an alternative to get_hgnc_name and is useful when
    HGNC name is not availabe (for instance, when the organism
    is not homo sapiens).

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
    protein_id = get_primary_id(protein_id)
    try:
        gene_name = um.uniprot_gene_name[protein_id]
        # There are cases when the entry is in the resource
        # table but the gene name is empty. Often this gene
        # name is actually available in the web service RDF
        # so here we return only if the gene name is not None
        # and not empty string.
        if gene_name:
            return gene_name
    except KeyError:
        pass
    if not web_fallback:
        return None

    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?gene a up:Gene .
            ?gene skos:prefLabel ?name .
            }
        """
    res = g.query(query)
    if res:
        gene_name = [r for r in res][0][0].toPython()
        if not gene_name:
            return None
        return gene_name
    return None


def get_gene_synonyms(protein_id):
    """Return a list of synonyms for the gene corresponding to a protein.

    Note that synonyms here also include the official gene name as
    returned by get_gene_name.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein to query

    Returns
    -------
    synonyms : list[str]
        The list of synonyms of the gene corresponding to the protein
    """
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?gene skos:altLabel | skos:prefLabel ?name .
            }
        """
    res = g.query(query)
    if res:
        return [r[0].toPython() for r in res]
    return None


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
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?name
        WHERE {
            ?gene :fullName | :shortName ?name .
            }
        """
    res = g.query(query)
    if res:
        return [r[0].toPython() for r in res]
    return None


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
    ret = []
    gene_syms = get_gene_synonyms(protein_id)
    if gene_syms:
        ret.extend(gene_syms)
    prot_syms = get_protein_synonyms(protein_id)
    if prot_syms:
        ret.extend(prot_syms)
    return ret


@lru_cache(maxsize=1000)
def get_sequence(protein_id):
    try:
        prim_ids = um.uniprot_sec[protein_id]
        protein_id = prim_ids[0]
    except KeyError:
        pass
    # Try to get the sequence from the downloaded sequence files
    seq = um.uniprot_sequences.get(protein_id)
    if seq is None:
        url = uniprot_url + '%s.fasta' % protein_id
        res = requests.get(url)
        res.raise_for_status()
        # res.text is Unicode
        lines = res.text.splitlines()
        seq = (''.join(lines[1:])).replace('\n','')
    return seq


def get_modifications(protein_id):
    g = query_protein(protein_id)
    if g is None:
        return None
    query = rdf_prefixes + """
        SELECT ?beg_pos ?comment
        WHERE {
            ?mod_res a up:Modified_Residue_Annotation .
            ?mod_res rdfs:comment ?comment .
            ?mod_res up:range ?range .
            ?range faldo:begin ?beg .
            ?range faldo:end ?end .
            ?beg a faldo:ExactPosition .
            ?beg faldo:position ?beg_pos .
            FILTER (?beg = ?end)
            }
        """
    res = g.query(query)
    mods = []
    for r in res:
        mod_pos = r[0].value
        # "Phosphothreonine; by autocatalysis"
        # "Phosphothreonine; by MAP2K1 and MAP2K2"
        # TODO: take into account the comment after the ;?
        mod_res = r[1].value.split(';')[0]
        mods.append((mod_res, mod_pos))
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
    return _is_organism(protein_id, 'RAT')


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


def get_mouse_id(human_protein_id):
    """Return the mouse UniProt ID given a human UniProt ID.

    Parameters
    ----------
    human_protein_id : str
        The UniProt ID of a human protein.

    Returns
    -------
    mouse_protein_id : str
        The UniProt ID of a mouse protein orthologous to the given human protein
    """
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
    return um.uniprot_length.get(protein_id)


@lru_cache(maxsize=10000)
def query_protein_xml(protein_id):
    """Retrieve the XML entry for a given protein.

    Some information is only available in the XML entry for UniProt
    proteins (not RDF), therefore this endpoint is necessary.

    Parameters
    ----------
    protein_id : str
        The UniProt ID of the protein to look up.

    Returns
    -------
    xml.etree.ElementTree
        An ElementTree representation of the XML entry for the
        protein.
    """
    try:
        prim_ids = um.uniprot_sec[protein_id]
        protein_id = prim_ids[0]
    except KeyError:
        pass
    url = uniprot_url + protein_id + '.xml'
    try:
        ret = requests.get(url)
        et = ElementTree.fromstring(ret.content)
    except Exception as e:
        return None
    return et


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
    et = query_protein_xml(protein_id)
    if et is None:
        return None
    function = et.find('up:entry/up:comment[@type="function"]/up:text',
                       namespaces={'up': 'http://uniprot.org/uniprot'})
    if function is None:
        return None
    return function.text


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
    tuple of int
        The beginning and end position of the signal peptide as a tuple
        of integers.
    """
    # Note, we use False here to differentiate from None
    entry = um.signal_peptide.get(protein_id, False)
    if entry is not False or not web_fallback:
        return entry
    et = query_protein_xml(protein_id)
    if et is None:
        return None, None
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
    return begin_pos, end_pos


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


class UniprotMapper(object):
    def __init__(self):
        self.initialized = False
        self.initialized_hmr = False
        self.initialized_seq = False
        self.initialized_refseq = False

    def initialize(self):
        maps = _build_uniprot_entries()
        (self._uniprot_gene_name, self._uniprot_mnemonic, \
         self._uniprot_mnemonic_reverse, self._uniprot_mgi,
         self._uniprot_rgd, self._uniprot_mgi_reverse,
         self._uniprot_rgd_reverse, self._uniprot_length,
         self._uniprot_reviewed, self._uniprot_signal_peptide) = maps

        self._uniprot_sec = _build_uniprot_sec()

        self.initialized = True

    def initialize_hmr(self):
        self._uniprot_human_mouse, self._uniprot_human_rat = \
            _build_human_mouse_rat()
        self.initialized_hmr = True

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
        if not self.initialized_hmr:
            self.initialize_hmr()
        return self._uniprot_human_mouse

    @property
    def uniprot_human_rat(self):
        if not self.initialized_hmr:
            self.initialize_hmr()
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
    def signal_peptide(self):
        if not self.initialized:
            self.initialize()
        return self._uniprot_signal_peptide


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
    uniprot_signal_peptide = {}
    uniprot_reviewed = set()
    with open(up_entries_file, 'r') as fh:
        csv_rows = csv.reader(fh, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        for row in csv_rows:
            up_id, gene_name, up_mnemonic, rgd, mgi, length, reviewed, \
                signal_peptide = row
            # Store the entry in the reviewed set
            if reviewed == 'reviewed':
                uniprot_reviewed.add(up_id)
            uniprot_gene_name[up_id] = gene_name
            uniprot_mnemonic[up_id] = up_mnemonic
            uniprot_mnemonic_reverse[up_mnemonic] = up_id
            uniprot_length[up_id] = int(length)
            if mgi:
                mgi_ids = mgi.split(';')
                if mgi_ids:
                    uniprot_mgi[up_id] = mgi_ids[0]
                    uniprot_mgi_reverse[mgi_ids[0]] = up_id
            if rgd:
                rgd_ids = rgd.split(';')
                if rgd_ids:
                    uniprot_rgd[up_id] = rgd_ids[0]
                    uniprot_rgd_reverse[rgd_ids[0]] = up_id
            uniprot_signal_peptide[up_id] = (None, None)
            if signal_peptide:
                match = re.match(r'SIGNAL (\d+) (\d+) ', signal_peptide)
                if match:
                    beg_pos, end_pos = match.groups()
                    uniprot_signal_peptide[up_id] = \
                        (int(beg_pos), int(end_pos))

    return (uniprot_gene_name, uniprot_mnemonic, uniprot_mnemonic_reverse,
            uniprot_mgi, uniprot_rgd, uniprot_mgi_reverse, uniprot_rgd_reverse,
            uniprot_length, uniprot_reviewed, uniprot_signal_peptide)


def _build_human_mouse_rat():
    hgnc_file = resource_manager.get_create_resource_file('hgnc')
    with open(hgnc_file, 'r') as fh:
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
    with open(hgnc_file, 'r') as fh:
        csv_rows = csv.reader(fh, delimiter='\t')
        # Skip the header row
        next(csv_rows)
        hgnc_ids = {}
        uniprot_ids = {}
        for row in csv_rows:
            hgnc_id = row[0][5:]
            hgnc_status = row[3]
            if hgnc_status == 'Approved':
                hgnc_name = row[1]
                hgnc_ids[hgnc_name] = hgnc_id
            # Uniprot
            uniprot_id = row[6]
            uniprot_ids[hgnc_id] = uniprot_id
    return hgnc_ids, uniprot_ids


def _build_uniprot_sec():
    # File containing secondary accession numbers mapped
    # to primary accession numbers
    sec_file = resource_manager.get_create_resource_file('upsec')
    uniprot_sec = {}
    lines = open(sec_file, 'rt').readlines()
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
    with open(refseq_uniprot_file, 'rt') as f:
        csvreader = csv.reader(f)
        for refseq_id, up_id in csvreader:
            if refseq_id not in refseq_up:
                refseq_up[refseq_id] = []
            refseq_up[refseq_id].append(up_id)
    return refseq_up


def load_fasta_sequences(seq_file, id_delimiter='|', id_index=1):
    sequences = {}
    with open(seq_file, 'rt') as f:
        lines = f.readlines()
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
