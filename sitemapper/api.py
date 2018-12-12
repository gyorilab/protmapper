from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from future.utils import python_2_unicode_compatible
import os
import pickle
import logging
import textwrap
from copy import deepcopy
from indra.statements import *
from indra.util import read_unicode_csv
from indra.config import has_config, get_config
from indra.databases import uniprot_client, hgnc_client, phosphosite_client

# Python 2
try:
    basestring
# Python 3
except:
    basestring = str

logger = logging.getLogger(__name__)


if has_config('SITEMAPPER_CACHE_PATH'):
    sitemapper_cache = get_config('SITEMAPPER_CACHE_PATH')
else:
    sitemapper_cache = None

class MappedStatement(object):
    """Information about a Statement found to have invalid sites.

    Parameters
    ----------
    original_stmt : :py:class:`indra.statements.Statement`
        The statement prior to mapping.
    mapped_mods : list of tuples
        A list of invalid sites, where each entry in the list has two
        elements: ((gene_name, residue, position), mapped_site).  If the
        invalid position was not found in the site map, mapped_site is
        None; otherwise it is a tuple consisting of (residue, position,
        comment). Note that some entries in the site map are curated *errors*,
        that is, sites that are known to be frequent misattributions
        to certain proteins. Such sites are mapped to tuples
        (None, None, comment).
    mapped_stmt : :py:class:`indra.statements.Statement`
        The statement after mapping. Note that if no information was found
        in the site map, it will be identical to the original statement.
    """
    def __init__(self, original_stmt, mapped_mods, mapped_stmt):
        self.original_stmt = original_stmt
        self.mapped_mods = mapped_mods
        self.mapped_stmt = mapped_stmt

    @python_2_unicode_compatible
    def __str__(self):
        if not self.mapped_mods:
            mm_str = str(self.mapped_mods)
        else:
            mm_ws = '\n' + (' ' * 17)
            mm_str = mm_ws.join([str(mm) for mm in self.mapped_mods])

        summary = textwrap.dedent("""
            MappedStatement:
                original_stmt: {0}
                mapped_mods: {1}
                mapped_stmt: {2}
            """)
        return summary.format(self.original_stmt, mm_str, self.mapped_stmt)

    def __repr__(self):
        return str(self)


class SiteMapper(object):
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
    def __init__(self, site_map=None, use_cache=False):
        self.site_map = site_map
        self.use_cache = use_cache
        self._cache = {}
        if self.use_cache:
            self._cache_path = sitemapper_cache
            if os.path.exists(self._cache_path):
                with open(self._cache_path, 'rb') as f:
                    self._cache = pickle.load(f)
                print("Loaded cache of length %d." % len(self._cache))
        self._sitecount = {}

    def __del__(self):
        try:
            if self.use_cache:
                import pickle
                with open(self._cache_path, 'wb') as f:
                    pickle.dump(self._cache, f, protocol=2)
        except:
            pass

    def map_sites(self, prot_id, prot_ns, site_list):
        return 1234


    def _check_agent_mod(self, agent, mods, do_methionine_offset=True,
                         do_orthology_mapping=True,
                         do_isoform_mapping=True):
        """Check an agent for invalid sites and look for mappings.

        Look up each modification site on the agent in Uniprot and then the
        site map.

        Parameters
        ----------
        agent : :py:class:`indra.statements.Agent`
            Agent to check for invalid modification sites.
        mods : list of :py:class:`indra.statements.ModCondition`
            Modifications to check for validity and map.
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
        list
            A list of invalid sites, where each entry in the list has two
            elements: ((gene_name, residue, position), mapped_site).  If the
            invalid position was not found in the site map, mapped_site is
            None; otherwise it is a tuple consisting of (residue, position,
            comment).
        """
        invalid_sites = []
        up_id = _get_uniprot_id(agent)
        # If the uniprot entry is not found, let it pass
        if not up_id:
            logger.debug("No uniprot ID for %s" % agent.name)
            return [] # Same effect as valid sites
        # Look up all of the modifications in uniprot, and add them to the list
        # of invalid sites if they are missing
        for old_mod in mods:
            # If no site information for this residue, skip
            if old_mod.position is None or old_mod.residue is None:
                continue
            site_key = (agent.name, old_mod.residue, old_mod.position)
            # Increase our count for this site
            self._sitecount[site_key] = self._sitecount.get(site_key, 0) + 1
            # First, check the cache to potentially avoid a costly sequence
            # lookup
            cached_site = self._cache.get(site_key)
            if cached_site is not None:
                if cached_site == 'VALID':
                    pass
                else:
                    invalid_sites.append((site_key, cached_site))
                continue
            # If not cached, continue
            # Look up the residue/position in uniprot
            site_valid = uniprot_client.verify_location(up_id,
                                                        old_mod.residue,
                                                        old_mod.position)
            # If it's not found in Uniprot, then look it up in the site map
            if site_valid:
                self._cache[site_key] = 'VALID'
                continue
            # Check the agent for a Uniprot ID
            up_id = agent.db_refs.get('UP')
            hgnc_id = agent.db_refs.get('HGNC')
            if not hgnc_id:
                logger.debug("No HGNC ID for %s, only curated sites will be "
                            "mapped" % agent.name)
            # NOTE: The following lookups can only be performed if the
            # Phosphosite Data is available.
            if phosphosite_client.has_data():
                # First, look for other entries in phosphosite for this protein
                # where this sequence position is legit (i.e., other isoforms)
                if do_isoform_mapping and up_id and hgnc_id:
                    human_pos = phosphosite_client.map_to_human_site(
                                  up_id, old_mod.residue, old_mod.position)
                    if human_pos:
                        mapped_site = (old_mod.residue, human_pos,
                                       'INFERRED_ALTERNATIVE_ISOFORM')
                        self._cache[site_key] = mapped_site
                        invalid_sites.append((site_key, mapped_site))
                        continue
                # Try looking for rat or mouse sites
                if do_orthology_mapping and up_id and hgnc_id:
                    # Get the mouse ID for this protein
                    up_mouse = uniprot_client.get_mouse_id(up_id)
                    # Get mouse sequence
                    human_pos = phosphosite_client.map_to_human_site(
                                  up_mouse, old_mod.residue, old_mod.position)
                    if human_pos:
                        mapped_site = (old_mod.residue, human_pos,
                                       'INFERRED_MOUSE_SITE')
                        self._cache[site_key] = mapped_site
                        invalid_sites.append((site_key, mapped_site))
                        continue
                    # Try the rat sequence
                    up_rat = uniprot_client.get_rat_id(up_id)
                    human_pos = phosphosite_client.map_to_human_site(
                                  up_rat, old_mod.residue, old_mod.position)
                    if human_pos:
                        mapped_site = (old_mod.residue, human_pos,
                                       'INFERRED_RAT_SITE')
                        self._cache[site_key] = mapped_site
                        invalid_sites.append((site_key, mapped_site))
                        continue
                # Check for methionine offset (off by one)
                if do_methionine_offset and up_id and hgnc_id:
                    try:
                        offset_pos = str(int(old_mod.position) + 1)
                    except ValueError:
                        logger.warning("Invalid position: %s" %
                                       old_mod.position)
                        continue
                    human_pos = phosphosite_client.map_to_human_site(
                                  up_id, old_mod.residue, offset_pos)
                    # If it's valid at the offset position, create the mapping
                    # and continue
                    if human_pos:
                        mapped_site = (old_mod.residue, human_pos,
                                       'INFERRED_METHIONINE_CLEAVAGE')
                        self._cache[site_key] = mapped_site
                        invalid_sites.append((site_key, mapped_site))
                        continue
            # Now check the site map
            mapped_site = self.site_map.get(site_key, None)
            # No entry in the site map; set site info to None
            if mapped_site is None:
                self._cache[site_key] = None
                invalid_sites.append((site_key, None))
            # Manually mapped in the site map
            else:
                self._cache[site_key] = mapped_site
                invalid_sites.append((site_key, mapped_site))
        return invalid_sites


def _get_uniprot_id(agent):
    """Get the Uniprot ID for an agent, looking up in HGNC if necessary.

    If the Uniprot ID is a list then return the first ID by default.
    """
    up_id = agent.db_refs.get('UP')
    hgnc_id = agent.db_refs.get('HGNC')
    if up_id is None:
        if hgnc_id is None:
            # If both UniProt and HGNC refs are missing we can't
            # sequence check and so don't report a failure.
            return None
        # Try to get UniProt ID from HGNC
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        # If this fails, again, we can't sequence check
        if up_id is None:
            return None
    # If the UniProt ID is a list then choose the first one.
    if not isinstance(up_id, basestring) and \
       isinstance(up_id[0], basestring):
        up_id = up_id[0]
    return up_id


def load_site_map(path):
    """Load the modification site map from a file.

    The site map file should be a comma-separated file with six columns::

        Gene: HGNC gene name
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
        A dict mapping tuples of the form `(gene, orig_res, orig_pos)` to a
        tuple of the form `(correct_res, correct_pos, comment)`, where `gene`
        is the string name of the gene (canonicalized to HGNC); `orig_res` and
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
        if not (row[0] and row[1] and row[2]):
            raise Exception("Entries in the key (gene, residue, position) "
                            "may not be empty.")
        correct_res = row[3].strip() if row[3] else None
        correct_pos = row[4].strip() if row[4] else None
        comment = row[5].strip() if row[5] else None
        site_map[(row[0].strip(), row[1].strip(), row[2].strip())] = \
                                (correct_res, correct_pos, comment)
    return site_map


default_site_map_path = os.path.join(os.path.dirname(__file__),
                             '../resources/curated_site_map.csv')

default_site_map = load_site_map(default_site_map_path)

default_mapper = SiteMapper(default_site_map)
"""A default instance of :py:class:`SiteMapper` that contains the site
information found in resources/curated_site_map.csv'."""

