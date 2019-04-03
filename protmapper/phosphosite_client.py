import csv
import logging
from os.path import dirname, abspath, join
from collections import namedtuple, defaultdict
from protmapper.resources import resource_manager


logger = logging.getLogger(__name__)


PhosphoSite = namedtuple('PhosphoSite',
                         ['GENE', 'PROTEIN', 'ACC_ID', 'HU_CHR_LOC', 'MOD_RSD',
                          'SITE_GRP_ID', 'ORGANISM', 'MW_kD', 'DOMAIN',
                          'SITE_7_AA', 'LT_LIT', 'MS_LIT', 'MS_CST', 'CST_CAT'])

PspMapping = namedtuple('PspMapping',
                        ['mapped_id', 'mapped_res', 'mapped_pos', 'motif',
                         'respos'])

_data_by_up = None
_data_by_site_grp = None
_has_data = None



# An explicit mapping of Uniprot IDs representing reference isoforms of their
# respective proteins in cases where there is ambiguity (i.e., other isoforms
# have IDs without a hyphen).
_iso_to_ref_map = {
    'J3KPC8': 'Q9Y2K2', # SIK3/QSK
    'Q5W9H1': 'Q14155', # ARHGEF7
    'A4QN18': 'O15027', # SEC16A
    'J3KNL6': 'O15027', # SEC16A
    'NP_001184222': 'Q16555', # DPYSL2/CRMP-2
}


# An easily indexed list of the reference isoforms mapped to, above
_ref_isoforms = set(_iso_to_ref_map.values())



def has_data():
    """Check if the PhosphoSite data is available and can be loaded.

    Returns
    -------
    bool
        True if the data can be loaded, False otherwise.
    """
    global _has_data
    if _has_data is None:
        try:
            _get_phospho_site_dataset()
            # If we succeeded without exception, then we set _has_data to True
            _has_data = True
        except Exception as e:
            logger.info("Could not load PhosphoSite data from file %s" %
                         resource_manager.get_resource_file('psp'))
            logger.info("Source Exception: %s" % e)
            _has_data = False
    return _has_data


def _get_phospho_site_dataset():
    """Read phosphosite data into dicts keyed by Uniprot ID and by site group.

    Returns
    -------
    tuple
        The first element of the tuple contains the PhosphoSite data keyed
        by Uniprot ID, the second element contains data keyed by site group.
        Both dicts have instances of the PhosphoSite namedtuple as values.
        If the PhosphoSite data file cannot be loaded, returns (None, None).
    """
    global _data_by_up
    global _data_by_site_grp
    phosphosite_data_file = resource_manager.get_create_resource_file('psp')
    if _data_by_up is None or _data_by_site_grp is None:
        with open(phosphosite_data_file, 'r') as fh:
            # Get the csv reader generator
            reader = csv.reader(fh, delimiter='\t')
            # Skip 4 rows
            for _ in range(4):
                next(reader)
            # Build up a dict by protein
            data_by_up = defaultdict(lambda: defaultdict(list))
            data_by_site_grp = defaultdict(list)
            for row in reader:
                site = PhosphoSite(*row)
                res_pos = site.MOD_RSD.split('-')[0]
                #res_pos = res_pos[1:] # DANGEROUS: lookup based on pos alone
                base_acc_id = site.ACC_ID.split('-')[0]
                data_by_up[site.ACC_ID][res_pos].append(site)
                # If the ID was isoform specific, add to the dict for the whole
                # protein
                if base_acc_id != site.ACC_ID:
                    data_by_up[base_acc_id][res_pos].append(site)
                # Catch the handful of isoforms that have a Uniprot ID without
                # the hyphen
                elif site.ACC_ID in _iso_to_ref_map:
                    ref_id = _iso_to_ref_map[site.ACC_ID]
                    data_by_up[ref_id][res_pos].append(site)
                # To catch additional cases, include an entry for the -1 base ID
                else:
                    data_by_up['%s-1' % base_acc_id] = data_by_up[base_acc_id]
                data_by_site_grp[site.SITE_GRP_ID].append(site)
            _data_by_up = data_by_up
            _data_by_site_grp = data_by_site_grp
    return (_data_by_up, _data_by_site_grp)


def map_to_human_site(up_id, mod_res, mod_pos):
    """Find site on human ref seq corresponding to (possibly non-human) site.

    Parameters
    ----------
    up_id : str
        Uniprot ID of the modified protein (generally human, rat, or mouse).
    mod_res : str
        Modified amino acid residue.
    mod_pos : str
        Amino acid sequence position.

    Returns
    -------
    str
        Returns amino acid position on the human reference sequence
        corresponding to the site on the given protein.
    """
    (data_by_up, data_by_site_grp) = _get_phospho_site_dataset()
    sites_for_up = data_by_up.get(up_id)
    # No info in Phosphosite for this Uniprot ID
    if not sites_for_up:
        return None
    site_info_list = sites_for_up.get('%s%s' % (mod_res, mod_pos))
    # DANGER: lookup based on pos alone
    # site_info_list = sites_for_up.get('%s' % mod_pos)
    # If this site doesn't exist for this protein, will return an empty list
    if not site_info_list:
        return None
    # At this point site_info_list contains a list of PhosphoSite objects
    # for the given protein with an entry for the given site, however, they
    # may be different isoforms; they may also not contain the reference
    # isoform; or they may only contain a single isoform.
    # If it's a single isoform, take it!
    if len(site_info_list) == 1:
        site_grp = site_info_list[0].SITE_GRP_ID
    # If there is more than one entry for this site, take the one from the
    # reference sequence, if it's in there (for example, a site that is present
    # in a mouse isoform as well as the mouse reference sequence)
    elif len(site_info_list) > 0:
        logger.debug('More than one entry in PhosphoSite for %s, site %s%s' %
                    (up_id, mod_res, mod_pos))
        ref_site_info = None
        si_by_grp = defaultdict(list)
        for si in site_info_list:
            logger.debug('\tSite info: acc_id %s, site_grp_id %s' %
                        (si.ACC_ID, si.SITE_GRP_ID))
            si_by_grp[si.SITE_GRP_ID].append(si)
            if si.ACC_ID == up_id:
                ref_site_info = si
        # If the reference sequence is not there but we have more than one
        # site_info in the list, this means there is more than one isoform
        # that has a phosphorylated site at the given position. The main thing
        # is whether there is more than one site group.
        if ref_site_info is None:
            # If there's only one unique site group, take it
            if len(si_by_grp) == 1:
                site_grp = site_info_list[0].SITE_GRP_ID
            # Otherwise, take the first
            else:
                logger.warning('More than one non-reference site group found '
                               'for (%s, %s, %s), not mapping' %
                               (up_id, mod_res, mod_pos))
                return None
        else:
            site_grp = ref_site_info.SITE_GRP_ID
    # Look up site group
    site_grp_list = data_by_site_grp[site_grp]
    # If an empty list, then return None (is unlikely to happen)
    if not site_grp_list:
        return None
    # Check for a human protein in the list
    human_sites = [s for s in site_grp_list if s.ORGANISM == 'human']
    if not human_sites:
        logger.debug("No corresponding human site for %s, choosing first: %s" %
                     (site_grp, site_grp_list[0].ORGANISM))
        #ref_site = site_grp_list[0]
        return None
    # If there are multiple isoforms, choose the base one
    # (no hyphen in the accession ID)
    elif len(human_sites) > 1:
        # In general we assume that multiple human sites only arise from
        # multiple isoforms of the same protein, which will share an accession
        # ID, and that only one of these will be the reference sequence (with
        # no hyphen).  However, in a small number of cases the isoform-specific
        # identifier will have a Uniprot ID without the hyphen, e.g. J3KPC8
        # for QSK iso5 (ref protein with ID Q9Y2K2).
        base_id_sites = [site for site in human_sites
                         if site.ACC_ID.find('-') == -1 and
                            site.ACC_ID not in _iso_to_ref_map]
        if base_id_sites:
            if len(base_id_sites) != 1:
                logger.warning("There is more than one apparent ref seq, "
                               "choosing first: %s" % base_id_sites)
            ref_site = base_id_sites[0]
        # There is no base ID site, i.e., all mapped sites are for specific
        # isoforms only, so skip it
        else:
            logger.info('Human isoform matches, but no human ref seq match '
                        'for %s %s %s; choosing first' %
                        (up_id, mod_res, mod_pos))
            ref_site = human_sites[0]
    # If there is only one human site, take it
    else:
        ref_site = human_sites[0]
    mapped_id = ref_site.ACC_ID
    human_site_str = ref_site.MOD_RSD.split('-')[0]
    human_res = human_site_str[0]
    human_pos = human_site_str[1:]
    motif, respos = _normalize_site_motif(ref_site.SITE_7_AA)
    pspmapping = PspMapping(mapped_id=mapped_id, mapped_res=human_res,
                            mapped_pos=human_pos,
                            motif=motif, respos=respos)
    return pspmapping


def sites_only(exclude_isoforms=False):
    """Return PhosphositePlus data as a flat list of proteins and sites.

    Parameters
    ----------
    exclude_isoforms : bool
        Whether to exclude sites for protein isoforms. Default is False
        (includes isoforms).

    Returns
    -------
    list of tuples
        Each tuple consists of (uniprot_id, residue, position).
    """
    sites = []
    (data_by_up, data_by_site_grp) = _get_phospho_site_dataset()
    for up_id, respos_dict in data_by_up.items():
        if exclude_isoforms and '-' in up_id:
            continue
        for respos in respos_dict.keys():
            res = respos[0]
            pos = respos[1:]
            sites.append((up_id, res, pos))
    return sites


def _normalize_site_motif(motif):
    """Normalize the PSP site motif to all caps with no underscores and return
    the preprocessed motif sequence and the position of the target residue
    in the motif (zero-indexed)."""
    no_underscores = motif.replace('_', '')
    offset = motif.find(no_underscores)
    respos = 7 - offset
    return (no_underscores.upper(), respos)
