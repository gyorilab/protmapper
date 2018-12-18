from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import os
import pickle
import logging
import textwrap
from copy import deepcopy
from indra.util import read_unicode_csv
from indra.databases import uniprot_client, hgnc_client
from protmapper import phosphosite_client
from protmapper.resources import resource_dir


# Python 2
try:
    basestring
# Python 3
except:
    basestring = str


logger = logging.getLogger(__name__)


valid_aas = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
             'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')


class MappedSite(object):
    """Represent details of a site that was mapped.

    Attributes
    ----------
    up_id : str
        The UniProt ID of the protein whose site was mapped.
    valid : bool
        True if the original site was valid with respect to the given
        protein, Falso otherwise.
    orig_res : str
        The original amino acid residue that was mapped.
    orig_pos : str
        The original amino acid position that was mapped.
    mapped_res : str
        The mapped amino acid residue.
    mapped_pos : str
        The mapped amino acid position.
    description : str
        A description of the mapping that was done, comes from a fixed
        set of codes of types of mapping that were performed.
    gene_name : str
        The standard (HGNC) gene name of the protein that was mapped.
    """
    def __init__(self, up_id, valid, orig_res, orig_pos, mapped_res=None,
                 mapped_pos=None, description=None, gene_name=None):
        self.up_id = up_id
        self.valid = valid
        self.orig_res = orig_res
        self.orig_pos = orig_pos
        self.mapped_res = mapped_res
        self.mapped_pos = mapped_pos
        self.description = description
        self.gene_name = gene_name

    def __repr__(self):
        return ("MappedSite(up_id='%s', valid=%s, orig_res='%s', "
                           "orig_pos='%s', mapped_res='%s', mapped_pos='%s', "
                           "description='%s', gene_name='%s')" %
                           (self.up_id, self.valid, self.orig_res,
                            self.orig_pos, self.mapped_res, self.mapped_pos,
                            self.description, self.gene_name))

    def __eq__(self, other):
        if (self.up_id == other.up_id and self.valid == other.valid and
            self.orig_res == other.orig_res and
            self.orig_pos == other.orig_pos and
            self.mapped_res == other.mapped_res and
            self.mapped_pos == other.mapped_pos and
            self.description == other.description and
            self.gene_name == other.gene_name):
            return True
        else:
            return False

    def __ne__(self, other):
        return not(self == other)

    def __hash__(self):
        return hash((self.up_id, self.orig_res, self.orig_pos, self.mapped_res,
                     self.mapped_pos, self.description, self.gene_name))

    def to_json(self):
        keys = ('up_id', 'valid', 'orig_res', 'orig_pos', 'mapped_res',
                'mapped_pos', 'description', 'gene_name')
        jd = {key: self.__dict__.get(key) for key in keys}
        return jd


class ProtMapper(object):
    """
    Use curated site information to standardize modification sites in stmts.

    Parameters
    ----------
    site_map : dict (as returned by :py:func:`load_site_map`)
        A dict mapping tuples of the form `(gene, orig_res, orig_pos)` to a
        tuple of the form `(correct_res, correct_pos, comment)`, where `gene`
        is the string name of the gene (canonicalized to HGNC); `orig_res` and
        `orig_pos` are the residue and position to be mapped; `correct_res` and
        `correct_pos` are the corrected residue and position, and `comment` is
        a string describing the reason for the mapping (species error, isoform
        error, wrong residue name, etc.).
    use_cache : Optional[bool]
        If True, the SITEMAPPER_CACHE_PATH from the config (or environment)
        is loaded and cached mappings are read and written to the given path.
        Otherwise, no cache is used. Default: False

    Examples
    --------
    Fixing site errors on both the modification state of an agent (MAP2K1) and
    the target of a Phosphorylation statement (MAPK1):

    >>> map2k1_phos = Agent('MAP2K1', db_refs={'UP':'Q02750'}, mods=[
    ... ModCondition('phosphorylation', 'S', '217'),
    ... ModCondition('phosphorylation', 'S', '221')])
    >>> mapk1 = Agent('MAPK1', db_refs={'UP':'P28482'})
    >>> stmt = Phosphorylation(map2k1_phos, mapk1, 'T','183')
    >>> (valid, mapped) = default_mapper.map_sites([stmt])
    >>> valid
    []
    >>> mapped  # doctest:+IGNORE_UNICODE
    [
    MappedStatement:
        original_stmt: Phosphorylation(MAP2K1(mods: (phosphorylation, S, 217), (phosphorylation, S, 221)), MAPK1(), T, 183)
        mapped_mods: (('MAP2K1', 'S', '217'), ('S', '218', 'off by one'))
                     (('MAP2K1', 'S', '221'), ('S', '222', 'off by one'))
                     (('MAPK1', 'T', '183'), ('T', '185', 'off by two; mouse sequence'))
        mapped_stmt: Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222)), MAPK1(), T, 185)
    ]
    >>> ms = mapped[0]
    >>> ms.original_stmt
    Phosphorylation(MAP2K1(mods: (phosphorylation, S, 217), (phosphorylation, S, 221)), MAPK1(), T, 183)
    >>> ms.mapped_mods # doctest:+IGNORE_UNICODE
    [(('MAP2K1', 'S', '217'), ('S', '218', 'off by one')), (('MAP2K1', 'S', '221'), ('S', '222', 'off by one')), (('MAPK1', 'T', '183'), ('T', '185', 'off by two; mouse sequence'))]
    >>> ms.mapped_stmt
    Phosphorylation(MAP2K1(mods: (phosphorylation, S, 218), (phosphorylation, S, 222)), MAPK1(), T, 185)
    """
    def __init__(self, site_map=None, use_cache=False, cache_path=None):
        if site_map is None:
            site_map = load_site_map(default_site_map_path)

        self.site_map = site_map
        # Set default cache path
        if cache_path is None:
            cache_path = os.path.join(resource_dir, 'sm_cache.pkl')
        self._cache_path = cache_path
        self.use_cache = use_cache
        self._cache = {}
        if self.use_cache:
            if os.path.exists(self._cache_path):
                with open(self._cache_path, 'rb') as f:
                    self._cache = pickle.load(f)
                print("Loaded cache of length %d from %s" %
                      (len(self._cache), self._cache_path))
            else:
                print("No cache found at %s, one will be created." %
                      self._cache_path)
        self._sitecount = {}

    def __del__(self):
        try:
            if self.use_cache:
                with open(self._cache_path, 'wb') as f:
                    pickle.dump(self._cache, f, protocol=2)
        except:
            pass

    def map_sitelist_to_human_ref(self, site_list, **kwargs):
        mapped_sites = []
        for ix, (prot_id, prot_ns, residue, position) in enumerate(site_list):
            logger.info("Mapping site %d of %d, cache size %d" %
                        (ix, len(site_list), len(self._cache)))
            try:
                ms = self.map_to_human_ref(prot_id, prot_ns, residue, position,
                                           **kwargs)
                mapped_sites.append(ms)
            except Exception as e:
                logger.error("Error occurred mapping site "
                             "(%s, %s, %s, %s): %s" %
                             (prot_id, prot_ns, residue, position, str(e)))
        return mapped_sites

    def map_to_human_ref(self, prot_id, prot_ns, residue, position,
                         do_methionine_offset=True,
                         do_orthology_mapping=True,
                         do_isoform_mapping=True):
        """Check an agent for invalid sites and look for mappings.

        Look up each modification site on the agent in Uniprot and then the
        site map.

        Parameters
        ----------
        prot_id : str
            A Uniprot ID or HGNC gene symbol for the protein.
        prot_ns : str
            One of 'uniprot' or 'hgnc' indicating the type of ID given.
        residue : str
            Residue to map on the protein to check for validity and map.
        position : str
            Position of the residue to check for validity and map.
        do_methionine_offset : boolean
            Whether to check for off-by-one errors in site position (possibly)
            attributable to site numbering from mature proteins after
            cleavage of the initial methionine. If True, checks the reference
            sequence for a known modification at 1 site position greater
            than the given one; if there exists such a site, creates the
            mapping. Default is True.
        do_orthology_mapping : boolean
            Whether to check sequence positions for known modification sites
            in mouse or rat sequences (based on PhosphoSitePlus data). If a
            mouse/rat site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.
        do_isoform_mapping : boolean
            Whether to check sequence positions for known modifications
            in other human isoforms of the protein (based on PhosphoSitePlus
            data). If a site is found that is linked to a site in the human
            reference sequence, a mapping is created. Default is True.

        Returns
        -------
        MappedSite
            The MappedSite object gives information on results of mapping the
            site. See :py:class:`protmapper.MappedSite` documentation for
            details.
        """
        # Check the protein ID and namespace
        if prot_id is None:
            raise ValueError("prot_id must not be None.")
        if prot_ns not in ('uniprot', 'hgnc'):
            raise ValueError("prot_ns must be either 'uniprot' or 'hgnc' (for "
                             "HGNC symbols)")
        # Make sure the sites are valid
        valid_res, valid_pos = _validate_site(residue, position)
        # Get Uniprot ID and gene name
        up_id = _get_uniprot_id(prot_id, prot_ns)
        gene_name = uniprot_client.get_gene_name(up_id)
        # If an HGNC ID was given and the uniprot entry is not found,
        # let it pass
        if up_id is None:
            assert prot_ns == 'hgnc' and prot_id is not None
            return MappedSite(None, True, residue, position,
                              description="NO_UNIPROT_ID",
                              gene_name=prot_id)
        site_key = (up_id, residue, position)
        # Increase our count for this site
        self._sitecount[site_key] = self._sitecount.get(site_key, 0) + 1
        # First, check the cache to potentially avoid a costly sequence
        # lookup
        cached_site = self._cache.get(site_key)
        if cached_site is not None:
            return cached_site
        # If not cached, continue
        # Look up the residue/position in uniprot
        site_valid = uniprot_client.verify_location(up_id, residue, position)
        # It's a valid site
        if site_valid:
            mapped_site = MappedSite(up_id, True, residue, position,
                                     mapped_res=residue, mapped_pos=position,
                                     description='VALID',
                                     gene_name=gene_name)
            self._cache[site_key] = mapped_site
            return mapped_site
        # NOTE: The following lookups can only be performed if the
        # Phosphosite Data is available.
        human_prot = uniprot_client.is_human(up_id)
        if phosphosite_client.has_data():
            # First, look for other entries in phosphosite for this protein
            # where this sequence position is legit (i.e., other isoforms)
            if do_isoform_mapping and up_id and human_prot:
                human_pos = phosphosite_client.map_to_human_site(
                              up_id, residue, position)
                if human_pos:
                    mapped_site = \
                         MappedSite(up_id, False, residue, position,
                                    mapped_res=residue, mapped_pos=human_pos,
                                    description='INFERRED_ALTERNATIVE_ISOFORM',
                                    gene_name=gene_name)
                    self._cache[site_key] = mapped_site
                    return mapped_site
            # Try looking for rat or mouse sites
            if do_orthology_mapping and up_id and human_prot:
                # Get the mouse ID for this protein
                up_mouse = uniprot_client.get_mouse_id(up_id)
                # Get mouse sequence
                human_pos = phosphosite_client.map_to_human_site(
                              up_mouse, residue, position)
                if human_pos:
                    mapped_site = \
                         MappedSite(up_id, False, residue, position,
                                    mapped_res=residue, mapped_pos=human_pos,
                                    description='INFERRED_MOUSE_SITE',
                                    gene_name=gene_name)
                    self._cache[site_key] = mapped_site
                    return mapped_site
                # Try the rat sequence
                up_rat = uniprot_client.get_rat_id(up_id)
                human_pos = phosphosite_client.map_to_human_site(
                              up_rat, residue, position)
                if human_pos:
                    mapped_site = (residue, human_pos, 'INFERRED_RAT_SITE')
                    mapped_site = \
                         MappedSite(up_id, False, residue, position,
                                    mapped_res=residue, mapped_pos=human_pos,
                                    description='INFERRED_RAT_SITE',
                                    gene_name=gene_name)
                    self._cache[site_key] = mapped_site
                    return mapped_site
            # Check for methionine offset (off by one)
            if do_methionine_offset and up_id and human_prot:
                offset_pos = str(int(position) + 1)
                human_pos = phosphosite_client.map_to_human_site(
                              up_id, residue, offset_pos)
                # If it's valid at the offset position, create the mapping
                # and continue
                if human_pos:
                    mapped_site = \
                         MappedSite(up_id, False, residue, position,
                                    mapped_res=residue, mapped_pos=human_pos,
                                    description='INFERRED_METHIONINE_CLEAVAGE',
                                    gene_name=gene_name)
                    self._cache[site_key] = mapped_site
                    return mapped_site
        # Now check the site map
        mapped_site = self.site_map.get(site_key, None)
        # No entry in the site map; set site info to None
        if mapped_site is None:
            mapped_site = MappedSite(up_id, False, residue, position,
                                     description='NO_MAPPING_FOUND',
                                     gene_name=gene_name)
            self._cache[site_key] = mapped_site
            return mapped_site
        # Manually mapped in the site map
        else:
            mapped_res, mapped_pos, description = mapped_site
            # Convert empty strings to None
            mapped_site = MappedSite(up_id, True, residue, position,
                                     mapped_res=mapped_res,
                                     mapped_pos=mapped_pos,
                                     description=description,
                                     gene_name=gene_name)
            self._cache[site_key] = mapped_site
            return mapped_site


def load_site_map(path):
    """Load the modification site map from a file.

    The site map file should be a comma-separated file with six columns::

        UniprotId: Uniprot ID of protein
        Gene: Gene name
        OrigRes: Original (incorrect) residue
        OrigPos: Original (incorrect) residue position
        CorrectRes: The correct residue for the modification
        CorrectPos: The correct residue position
        Comment: Description of the reason for the error.

    Parameters
    ----------
    path : string
        Path to the tab-separated site map file.

    Returns
    -------
    dict
        A dict mapping tuples of the form `(uniprot_id, orig_res, orig_pos)` to
        a tuple of the form `(correct_res, correct_pos, comment)`, where
        `uniprot_id` is the Uniprot ID of the protein; `orig_res` and
        `orig_pos` are the residue and position to be mapped; `correct_res` and
        `correct_pos` are the corrected residue and position, and `comment` is
        a string describing the reason for the mapping (species error, isoform
        error, wrong residue name, etc.).
    """
    site_map = {}
    maprows = read_unicode_csv(path)
    # Skip the header line
    next(maprows)
    for row in maprows:
        # Don't allow empty entries in the key section
        if not (row[0] and row[2] and row[3]):
            raise Exception("Entries in the key (gene, residue, position) "
                            "may not be empty.")
        correct_res = row[4].strip() if row[4] else None
        correct_pos = row[5].strip() if row[5] else None
        comment = row[6].strip() if row[6] else None
        site_map[(row[0].strip(), row[2].strip(), row[3].strip())] = \
                                (correct_res, correct_pos, comment)
    return site_map


def _validate_site(residue, position):
    if residue is None:
        raise ValueError('residue cannot be None')
    if position is None:
        raise ValueError('position cannot be None')
    # Check that the residue is a valid amino acid
    if residue not in valid_aas:
        raise ValueError('Residue %s not a valid amino acid' % residue)
    # Next make sure that the position is a valid position
    try:
        pos_int = int(position)
        pos_str = str(pos_int)
    except ValueError:
        raise ValueError('Position %s not a valid sequence position.'
                         % position)
    # Site appears valid, make a Site object
    return (residue, pos_str)


def _get_uniprot_id(prot_id, prot_ns):
    """Get the Uniprot ID for an agent, looking up in HGNC if necessary.

    If the Uniprot ID is a list then return the first ID by default.
    """
    # If the ID is a Uniprot ID, then we're done
    if prot_ns == 'uniprot':
        up_id = prot_id
    # Otherwise, we get the Uniprot ID from the HGNC name
    elif prot_ns == 'hgnc':
        # Get the HGNC ID
        hgnc_id = hgnc_client.get_hgnc_id(prot_id)
        if hgnc_id is None:
            return None
        # Try to get UniProt ID from HGNC
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        # If the UniProt ID is a list then choose the first one.
        if up_id is not None and \
           (not isinstance(up_id, basestring) and
            isinstance(up_id[0], basestring)):
            up_id = up_id[0]
    return up_id


default_site_map_path = os.path.join(os.path.dirname(__file__),
                             '../resources/curated_site_map.csv')

default_site_map = load_site_map(default_site_map_path)

default_mapper = ProtMapper(default_site_map)
"""A default instance of :py:class:`ProtMapper` that contains the site
information found in resources/curated_site_map.csv'."""

