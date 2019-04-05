import os
import csv
import pickle
import logging
from requests.exceptions import HTTPError
from protmapper.resources import resource_dir
from protmapper import phosphosite_client, uniprot_client


# Python 2
try:
    basestring
# Python 3
except:
    basestring = str


logger = logging.getLogger(__name__)


valid_aas = ('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
             'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')


class InvalidSiteException(Exception):
    pass


class MappedSite(object):
    """Represent details of a site that was mapped.

    Attributes
    ----------
    up_id : str
        The UniProt ID of the protein whose site was mapped.
    error_code : str or None
        One of several strings indicating an error in retrieving the protein
        sequence, or None if there was no error. Error codes include
        'NO_UNIPROT_ID' if the given gene name could not be converted into
        a Uniprot ID; 'UNIPROT_HTTP_NOT_FOUND' if the given Uniprot ID resulted
        in a 404 Not Found error from the Uniprot web service; or
        'UNIPROT_HTTP_OTHER' if it was any other type of Uniprot web service
        error. If the error code is not None, the `orig_res` and `orig_pos`
        fields will be set (based on the query arguments) but all other fields
        will be None.
    valid : bool
        True if the original site was valid with respect to the given
        protein, False otherwise. Further, in case of an error (if error_code is
        not None), it is set to None.
    orig_res : str
        The original amino acid residue that was mapped.
    orig_pos : str
        The original amino acid position that was mapped.
    mapped_id : str
        The Uniprot ID for the protein containing the mapped site. If `up_id`
        is the Uniprot ID for the human reference sequence, in most cases
        this will match; however, exceptions will occur if the site position
        refers to a site that is unique to a particular isoform.
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
    def __init__(self, up_id, valid, orig_res, orig_pos,
                 error_code=None, mapped_id=None, mapped_res=None,
                 mapped_pos=None, description=None, gene_name=None):
        self.up_id = up_id
        self.error_code = error_code
        self.valid = valid
        self.orig_res = orig_res
        self.orig_pos = orig_pos
        self.mapped_id = mapped_id
        self.mapped_res = mapped_res
        self.mapped_pos = mapped_pos
        self.description = description
        self.gene_name = gene_name

    def __repr__(self):
        quote_args = lambda args: tuple([a if a in (None, True, False)
                                           else ("'%s'" % a) for a in args])
        return ("MappedSite(up_id=%s, error_code=%s, valid=%s, "
                    "orig_res=%s, orig_pos=%s, mapped_id=%s, mapped_res=%s, "
                    "mapped_pos=%s, description=%s, gene_name=%s)" %
                quote_args([self.up_id, self.error_code, self.valid,
                            self.orig_res, self.orig_pos, self.mapped_id,
                            self.mapped_res, self.mapped_pos, self.description,
                            self.gene_name]))

    def __eq__(self, other):
        if (self.up_id == other.up_id and self.valid == other.valid and
            self.error_code == other.error_code and
            self.orig_res == other.orig_res and
            self.orig_pos == other.orig_pos and
            self.mapped_id == other.mapped_id and
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
        return hash((self.up_id, self.error_code, self.valid, self.orig_res,
                     self.orig_pos, self.mapped_id, self.mapped_res,
                     self.mapped_pos, self.description, self.gene_name))

    def to_json(self):
        keys = ('up_id', 'error_code', 'valid', 'orig_res', 'orig_pos',
                'mapped_id', 'mapped_res', 'mapped_pos', 'description',
                'gene_name')
        jd = {key: self.__dict__.get(key) for key in keys}
        return jd

    def not_invalid(self):
        """Return True if the original site is not known to be invalid.

        Returns
        -------
        bool
            True if the original site is valid or if there is an error code,
            which implicitly means that the validity of the original site could
            not be established. False otherwise.
        """
        return self.valid or self.error_code is not None

    def has_mapping(self):
        """Return True if the original site was mapped successfully.

        Returns
        -------
        bool
            True if a mapping was successfully obtained for the site, False
            otherwise.
        """
        return (not self.not_invalid()) and \
            (self.mapped_pos is not None and self.mapped_res is not None)


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

    def save_cache(self):
        with open(self._cache_path, 'wb') as f:
            pickle.dump(self._cache, f)

    def __del__(self):
        try:
            if self.use_cache:
                self.save_cache()
        except:
            pass

    def map_sitelist_to_human_ref(self, site_list, **kwargs):
        """Return a list of mapped sites for a list of input sites.

        Parameters
        ----------
        site_list : list of tuple
            Each tuple in the list consists of the following entries:
            (prot_id, prot_ns, residue, position).

        Returns
        -------
        list of :py:class:`protmapper.api.MappedSite`
            A list of MappedSite objects, one corresponding to each site in
            the input list.
        """
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
            site. See :py:class:`protmapper.api.MappedSite` documentation for
            details.
        """
        # Check the protein ID and namespace
        if prot_id is None:
            raise ValueError("prot_id must not be None.")
        if prot_ns not in ('uniprot', 'hgnc'):
            raise ValueError("prot_ns must be either 'uniprot' or 'hgnc' (for "
                             "HGNC symbols)")
        # Get Uniprot ID and gene name
        up_id = _get_uniprot_id(prot_id, prot_ns)
        # If an HGNC ID was given and the uniprot entry is not found, flag
        # as error
        if up_id is None:
            assert prot_ns == 'hgnc' and prot_id is not None
            return MappedSite(None, None, residue, position,
                              gene_name=prot_id, error_code='NO_UNIPROT_ID')
        # Make sure the sites are proper amino acids/positions
        try:
            valid_res, valid_pos = _validate_site(residue, position)
        except InvalidSiteException as ex:
            return MappedSite(up_id, None, residue, position,
                              error_code='INVALID_SITE',
                              description=str(ex))
        # Get the gene name from Uniprot
        gene_name = uniprot_client.get_gene_name(up_id)
        site_key = (up_id, residue, position)
        # First, check the cache to potentially avoid a costly sequence
        # lookup
        cached_site = self._cache.get(site_key)
        if cached_site is not None:
            return cached_site
        # If not cached, continue
        # Look up the residue/position in uniprot
        try:
            site_valid = uniprot_client.verify_location(up_id, residue,
                                                        position)
        except HTTPError as ex:
            if ex.response.status_code == 404:
                error_code = 'UNIPROT_HTTP_NOT_FOUND'
            else:
                error_code = 'UNIPROT_HTTP_OTHER'
            # Set error_code; valid will set to None, not True/False
            mapped_site = MappedSite(up_id, None, residue, position,
                                     error_code=error_code)
            return mapped_site
        # It's a valid site
        if site_valid:
            mapped_site = MappedSite(up_id, True, residue, position,
                                     description='VALID',
                                     gene_name=gene_name)
            self._cache[site_key] = mapped_site
            return mapped_site
        # If it's not a valid site, check the site map first
        curated_site = self.site_map.get(site_key, None)
        # Manually mapped in the site map
        if curated_site is not None:
            mapped_res, mapped_pos, description = curated_site
            mapped_site = MappedSite(up_id, False, residue, position,
                                     mapped_id=up_id,
                                     mapped_res=mapped_res,
                                     mapped_pos=mapped_pos,
                                     description=description,
                                     gene_name=gene_name)
            self._cache[site_key] = mapped_site
            return mapped_site

        # There is no manual mapping, next we try to see if UniProt
        # reports a signal peptide that could be responsible for the position
        # being shifted
        signal_peptide = uniprot_client.get_signal_peptide(up_id, False)
        # If there is valid signal peptide information from UniProt
        if signal_peptide and signal_peptide[0] == 1 and \
                signal_peptide[1] is not None:
            offset_pos = str(int(position) + signal_peptide[1])
            # Check to see if the offset position is known to be phosphorylated
            mapped_site = self.get_psp_mapping(
                                up_id, up_id, gene_name, residue, position,
                                offset_pos, 'SIGNAL_PEPTIDE_REMOVED')
            if mapped_site:
                return mapped_site
        # ...there's no manually curated site or signal peptide, so do mapping
        # via PhosphoSite if the data is available:
        human_prot = uniprot_client.is_human(up_id)
        if phosphosite_client.has_data():
            # First, look for other entries in phosphosite for this protein
            # where this sequence position is legit (i.e., other isoforms)
            if do_isoform_mapping and up_id and human_prot:
                mapped_site = self.get_psp_mapping(
                        up_id, up_id, gene_name, residue, position, position,
                        'INFERRED_ALTERNATIVE_ISOFORM')
                if mapped_site:
                    return mapped_site
            # Try looking for rat or mouse sites
            if do_orthology_mapping and up_id and human_prot:
                # Get the mouse ID for this protein
                up_mouse = uniprot_client.get_mouse_id(up_id)
                # Get mouse sequence
                mapped_site = self.get_psp_mapping(
                                    up_id, up_mouse, gene_name, residue,
                                    position, position, 'INFERRED_MOUSE_SITE')
                if mapped_site:
                    return mapped_site
                # Try the rat sequence
                up_rat = uniprot_client.get_rat_id(up_id)
                mapped_site = self.get_psp_mapping(
                                    up_id, up_rat, gene_name, residue, position,
                                    position, 'INFERRED_RAT_SITE')
                if mapped_site:
                    return mapped_site
            # Check for methionine offset (off by one)
            if do_methionine_offset and up_id and human_prot:
                offset_pos = str(int(position) + 1)
                mapped_site = self.get_psp_mapping(
                                    up_id, up_id, gene_name, residue, position,
                                    offset_pos, 'INFERRED_METHIONINE_CLEAVAGE')
                if mapped_site:
                    return mapped_site
        # If we've gotten here, the entry is 1) not in the site map, and
        # 2) we either don't have PSP data or no mapping was found using PSP
        mapped_site = MappedSite(up_id, False, residue, position,
                                 description='NO_MAPPING_FOUND',
                                 gene_name=gene_name)
        self._cache[site_key] = mapped_site
        return mapped_site

    def get_psp_mapping(self, orig_id, query_id, gene_name, res, pos,
                        query_pos, mapping_code):
        """
        Wrapper around Phosphosite queries that performs peptide remapping.

        The function is called with a uniprot ID, residue, and position
        combination that is used to query the phosphosite_client for a valid
        corresponding site on the human reference protein. The `mapping_code`
        is provided by the caller to indicate the type of mapping being
        attempted (e.g., human isoform, mouse, rat, methionine). If a valid
        mapping is obtained, this is the error code that is applied.  If a
        valid mapping is obtained but it is for a human isoform, this indicates
        that the queried site exists only on a human isoform and not on the
        human reference protein, and the code `ISOFORM_SPECIFIC_SITE` is used.
        If the site returned by the phosphosite_client is at a position that
        does not match the Uniprot reference sequence (which can happen when
        the queried site and the PhosphositePlus protein sequences both exclude
        the initial methionine), the site is remapped to the Uniprot reference
        sequence using the peptide information for the site in PhosphositePlus.
        In these cases, the mapping code `REMAPPED_FROM_PSP_SEQUENCE` is used.

        Parameters
        ----------
        orig_id : str
            Original Uniprot ID of the protein to be mapped.
        query_id : str
            Uniprot ID of the protein being queried for sites. This may differ
            from `orig_id` if the orthologous mouse or rat protein is being
            checked for sites.
        gene_name : str
            Gene name of the protein.
        res : str
            Residue of the site to be mapped.
        pos : str
            Position of the site to be mapped.
        query_pos : str
            Position being queried for a mapping. This differs from `pos`
            when off-by-one (methionine) errors are being checked.
        mapping_code : str
            Mapping code to apply in case of a successful mapping, e.g.
            `INFERRED_ALTERNATIVE_ISOFORM`, `INFERRED_MOUSE_SITE`, etc.

        Returns
        -------
        MappedSite or None
            MappedSite object containing the mapping, or None indicating
            that no mapping was found.
        """
        pspmapping = phosphosite_client.map_to_human_site(query_id, res,
                                                          query_pos)
        # If no mapping, return None
        if pspmapping is None:
            return None
        # If there is a mapping, check to make sure that it is valid wrt to the
        # reference sequence
        human_pos = pspmapping.mapped_pos
        # Check if the site mapped from PSP is valid in the Uniprot sequence
        # for the ID that we're interested in
        site_valid = uniprot_client.verify_location(pspmapping.mapped_id,
                                  pspmapping.mapped_res, pspmapping.mapped_pos)
        # If the mapped site is valid, we're done!
        if site_valid:
            # If the residue is different, change the code accordingly
            mapped_site = MappedSite(orig_id, False, res, pos,
                              mapped_id=pspmapping.mapped_id,
                              mapped_res=pspmapping.mapped_res,
                              mapped_pos=human_pos,
                              description=mapping_code, gene_name=gene_name)
        else:
            # If mapped site is invalid, attempt to re-map based on the seq
            updated_pos = ProtMapper.map_peptide(orig_id, pspmapping.motif,
                                                 pspmapping.respos)
            # If the re-mapping fails, we give up
            if updated_pos is None:
                return None
            # Otherwise, we update to the mapped position
            updated_pos_1x = str(updated_pos + 1)
            mapped_site = MappedSite(orig_id, False, res, pos,
                              mapped_id=pspmapping.mapped_id,
                              mapped_res=pspmapping.mapped_res,
                              mapped_pos=updated_pos_1x, # Switch to 1-indexed
                              description='REMAPPED_FROM_PSP_SEQUENCE',
                              gene_name=gene_name)
        site_key = (orig_id, res, pos)
        self._cache[site_key] = mapped_site
        return mapped_site

    @staticmethod
    def motif_from_position(up_id, pos, window=7):
        seq = uniprot_client.get_sequence(up_id)
        return ProtMapper.motif_from_position_seq(seq, pos, window)

    @staticmethod
    def motif_from_position_seq(seq, pos, window=7):
        # Validate that the position is an integer
        pos_str = str(pos)
        pos_int = int(pos)
        end_ix = pos_int + window if pos_int + window < len(seq) else len(seq)
        start_ix = pos_int - window - 1 if pos_int - window > 0 else 0
        motif = seq[start_ix:end_ix]
        # Get the residue position, which will be the same as the window
        # size unless the start_ix is something other than the pos - window
        respos = pos_int - start_ix
        return (motif, respos)

    @staticmethod
    def map_peptide(target_up_id, peptide, pos):
        seq = uniprot_client.get_sequence(target_up_id)
        peptide_start = seq.find(peptide)
        # Peptide not found in sequence, return None
        if peptide_start == -1:
            return None
        target_pos = peptide_start + pos
        return target_pos

    @staticmethod
    def map_peptide_to_human_ref(prot_id, prot_ns, peptide, site_pos):
        """Return a mapped site for a given peptide.

        Parameters
        ----------
        prot_id : str
            A Uniprot ID or HGNC gene symbol for the protein.
        prot_ns : str
            One of 'uniprot' or 'hgnc' indicating the type of ID given.
        peptide : str
            A string of amino acid symbols representing a peptide.
        site_pos : str
            A site position within the peptide. Note: site_pos is 1-indexed.

        Returns
        -------
        MappedSite
            The MappedSite object gives information on results of mapping the
            site. See :py:class:`protmapper.api.MappedSite` documentation for
            details.
        """
        # Get the uniprot ID for the gene
        # Check the protein ID and namespace
        if prot_id is None:
            raise ValueError("prot_id must not be None.")
        if prot_ns not in ('uniprot', 'hgnc'):
            raise ValueError("prot_ns must be either 'uniprot' or 'hgnc' (for "
                             "HGNC symbols)")
        if prot_ns  == 'uniprot' and len(prot_id.split('-')) != 1 and \
                                                prot_id.split('-')[1] != '1':
            raise ValueError("Protein ID passed in appears to be a "
                             "non-reference isoform: %s" % prot_id)
        # Get Uniprot ID and gene name
        up_id = _get_uniprot_id(prot_id, prot_ns)
        # If an HGNC ID was given and the uniprot entry is not found, flag
        # as error
        if up_id is None:
            assert prot_ns == 'hgnc' and prot_id is not None
            return MappedSite(None, None, None, None,
                              gene_name=prot_id, error_code='NO_UNIPROT_ID')
        # Get the gene name from Uniprot
        gene_name = uniprot_client.get_gene_name(up_id)
        mapped_pos = ProtMapper.map_peptide(up_id, peptide, site_pos)
        ms = MappedSite(up_id=up_id, valid=None, orig_res=None, orig_pos=None,
                        error_code=None, description=None, gene_name=gene_name)
        if mapped_pos is None:
            ms.valid = False
        else:
            ms.valid = True
            ms.mapped_id = up_id
            ms.mapped_res = peptide[site_pos - 1]
            ms.mapped_pos = str(mapped_pos)
        return ms


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
    with open(path, 'r') as fh:
        maprows = csv.reader(fh)
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
        raise InvalidSiteException('residue cannot be None')
    if position is None:
        raise InvalidSiteException('position cannot be None')
    # Check that the residue is a valid amino acid
    if residue not in valid_aas:
        raise InvalidSiteException('Residue %s not a valid amino acid' %
                                   residue)
    # Next make sure that the position is a valid position
    try:
        pos_int = int(position)
        pos_str = str(pos_int)
    # Catch the ValueError and re-raise as an InvalidSiteException
    except ValueError:
        raise InvalidSiteException('Position %s not a valid sequence position.'
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
        hgnc_id = hgnc_ids.get(prot_id)
        if not hgnc_id:
            return None
        # Try to get UniProt ID from HGNC
        up_id = uniprot_ids.get(hgnc_id)
        # If the UniProt ID is a list then choose the first one.
        if up_id and (not isinstance(up_id, basestring) and
                      isinstance(up_id[0], basestring)):
            up_id = up_id[0]
    return up_id


default_site_map_path = os.path.join(os.path.dirname(__file__),
                                     'curated_site_map.csv')

default_site_map = load_site_map(default_site_map_path)

default_mapper = ProtMapper(default_site_map)
"""A default instance of :py:class:`ProtMapper` that contains the site
information found in resources/curated_site_map.csv'."""


hgnc_ids, uniprot_ids = uniprot_client._build_hgnc_mappings()
