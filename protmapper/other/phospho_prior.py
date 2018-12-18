import random
import numpy as np
from indra.tools import assemble_corpus as ac
from itertools import chain
import pickle
from collections import defaultdict
from indra.util import read_unicode_csv, write_unicode_csv
from matplotlib_venn import venn2, venn3
from matplotlib import pyplot as plt
from indra.databases import hgnc_client, uniprot_client
from indra.statements import *
from indra.db.query_db_stmts import by_gene_role_type
import synapseclient
from indra.databases import omnipath_client
from collections import Counter

def get_ids(hgnc_name):
    hgnc_id = hgnc_client.get_hgnc_id(hgnc_name)
    up_id = hgnc_client.get_uniprot_id(hgnc_id)
    return {'HGNC': hgnc_id, 'UP': up_id}


def load_ov_sites(use_mapped=True):
    ov_sites = set([])
    mismatch_count = 0
    for row in read_unicode_csv('mapped_peptides.txt', delimiter='\t',
                                skiprows=1):
        # Check for uniprot ID
        gene_name = row[2]
        up_id = row[3]
        orig_site = row[4]
        mapped_site = row[-1]
        if orig_site != mapped_site:
            mismatch_count += 1
        if not up_id:
            continue
        if use_mapped:
            ov_sites.add((gene_name, mapped_site))
        else:
            ov_sites.add((gene_name, orig_site))
    print("%d mismatched sites of %d" % (mismatch_count, len(ov_sites)))
    return ov_sites


def load_prior_from_synapse(synapse_id='syn11273504', peptide_specific=True):
    syn = synapseclient.Synapse()
    syn.login()
    # Obtain a pointer and download the data
    syn_data = syn.get(synapse_id)
    prior = {}
    for row in read_unicode_csv(syn_data.path, delimiter='\t'):
        sub_name, site_info = row[0].split(':')
        res = site_info[0]
        pos = site_info[1:]
        gene_list = row[1].split(',')
        if peptide_specific:
            prior[site_info] = gene_list
        else:
            prior[sub_name] = gene_list
    return prior


def load_statements_from_synapse(synapse_id='syn11273504'):
    syn = synapseclient.Synapse()
    syn.login()
    # Obtain a pointer and download the data
    syn_data = syn.get(synapse_id)
    stmts = []
    for row in read_unicode_csv(syn_data.path, delimiter='\t'):
        sub_name, site_info = row[0].split(':')
        res = site_info[0]
        pos = site_info[1:]
        gene_list = row[1].split(',')
        for enz_name in gene_list:
            enz = Agent(enz_name, db_refs=get_ids(enz_name))
            sub = Agent(sub_name, db_refs=get_ids(sub_name))
            stmt = Phosphorylation(enz, sub, res, pos)
            stmts.append(stmt)
    stmts = ac.map_sequence(stmts)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_genes_only(stmts, specific_only=True)
    return stmts


def get_omnipath_stmts():
    stmts = omnipath_client.get_all_modifications()
    phos_stmts = ac.filter_by_type(stmts, Phosphorylation)
    dephos_stmts = ac.filter_by_type(stmts, Dephosphorylation)
    stmts = phos_stmts + dephos_stmts
    stmts = ac.map_sequence(stmts)
    stmts = ac.filter_human_only(stmts)
    #stmts = ac.filter_genes_only(stmts, specific_only=True)
    return stmts


def get_indra_phos_stmts():
    stmts = by_gene_role_type(stmt_type='Phosphorylation')
    stmts += by_gene_role_type(stmt_type='Dephosphorylation')
    stmts = ac.map_grounding(stmts)
    # Expand families before site mapping
    stmts = ac.expand_families(stmts)
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.map_sequence(stmts)
    ac.dump_statements(stmts, 'sources/indra_phos_sitemap.pkl')
    stmts = ac.run_preassembly(stmts, poolsize=4,
                               save='sources/indra_phos_stmts_pre.pkl')
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_genes_only(stmts, specific_only=True)
    ac.dump_statements(stmts, 'sources/indra_phos_stmts.pkl')
    return stmts


def get_indra_reg_act_stmts():
    try:
        stmts = ac.load_statements('sources/indra_reg_act_stmts.pkl')
        return stmts
    except:
        pass
    stmts = []
    for stmt_type in ('Activation', 'Inhibition', 'ActiveForm'):
        print("Getting %s statements from INDRA DB" % stmt_type)
        stmts += by_gene_role_type(stmt_type=stmt_type)
    stmts = ac.map_grounding(stmts, save='sources/indra_reg_act_gmap.pkl')
    stmts = ac.filter_grounded_only(stmts)
    stmts = ac.run_preassembly(stmts, poolsize=4,
                               save='sources/indra_reg_act_pre.pkl')
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_genes_only(stmts, specific_only=True)
    ac.dump_statements(stmts, 'sources/indra_reg_act_stmts.pkl')
    return stmts


def get_pc_reg_act_stmts():
    pass


def add_regulators(reg_stmts, prior, max_features=100):
    # Build a dict of regulators for each gene
    reg_dict = {}
    for stmt_ix, stmt in enumerate(reg_stmts):
        if not stmt.subj:
            continue
        subj = stmt.subj.name
        obj = stmt.obj.name
        if obj not in reg_dict:
            reg_dict[obj] = set([subj])
        else:
            reg_dict[obj].add(subj)

    ext_prior = {}
    for site, gene_set in prior.items():
        if len(gene_set) >= max_features:
            ext_prior[site] = set(list(gene_set)[0:max_features])

        reg_gene_set = set()
        for gene in gene_set:
            regs = reg_dict.get(gene)
            if regs:
                reg_gene_set |= regs
        if not regs:
            ext_prior[site] = gene_set
        else:
            cur_prior_length = len(gene_set)
            num_remaining_features = max_features - cur_prior_length
            shuffled_regs = list(regs)
            random.shuffle(shuffled_regs)
            added_features = set(shuffled_regs[0:num_remaining_features])
            ext_prior[site] = gene_set | added_features
    return ext_prior


def get_phosphosite_stmts():
    stmts = ac.load_statements('sources/phosphosite_stmts.pkl')
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_genes_only(stmts, specific_only=True)
    return stmts


def save_indra_db_stmts(stmts):
    csv_rows = [('KINASE', 'KINASE_TEXT', 'SUBSTRATE', 'SUBSTRATE_TEXT',
                 'RESIDUE', 'POSITION', 'SOURCE', 'DIRECT', 'PMID', 'SENTENCE')]
    for s in stmts:
        for e in s.evidence:
            is_direct = 'True' if e.epistemics.get('direct') else 'False'
            csv_rows.append((s.enz.name, s.enz.db_refs.get('TEXT'),
                             s.sub.name, s.sub.db_refs.get('TEXT'),
                             s.residue, s.position, e.source_api,
                             is_direct, e.pmid, e.text))
    write_unicode_csv('indra_phosphosites.csv', csv_rows)


def plot_overlap(indra_sites, ps_sites, nk_sites):
    print('INDRA: %d unique sites' % len(indra_sites))
    print('Phosphosite: %d unique sites' % len(phos_sites))
    print('Intersection: %d sites' % len(indra_sites.intersection(phos_sites)))
    print('Sites in INDRA (not Phosphosite): %d sites' %
            len(indra_sites.difference(phos_sites)))
    print('Sites in Phosphosite (not INDRA): %d sites' %
            len(phos_sites.difference(indra_sites)))

    indra_only = indra_sites.difference(phos_sites)

    plt.ion()
    plt.figure()
    venn3((indra_sites, phos_sites, nk_sites),
              set_labels=('REACH/INDRA', 'PhosphoSite', 'NetworKIN'))
    plt.savefig('kinase_substrate_overlap.pdf')


def to_prior(stmts):
    prior = {}
    for stmt in stmts:
        key = '%s:%s%s' % (stmt.sub.name, stmt.residue, stmt.position)
        if key not in prior:
            prior[key] = set([stmt.enz.name])
        else:
            prior[key].add(stmt.enz.name)
    return prior


def to_nonspec_prior(stmts):
    prior = {}
    for stmt in stmts:
        if stmt.enz is None or stmt.sub is None:
            continue
        enz = stmt.enz.name
        sub = stmt.sub.name
        if sub not in prior:
            prior[sub] = set([enz])
        else:
            prior[sub].add(enz)
    return prior


def save_prior(prior, prior_name):
    with open(prior_name, 'wt') as f:
        for gene_key in sorted(prior.keys()):
            enzyme_list = ','.join(prior[gene_key])
            f.write('%s\t%s\n' % (gene_key, enzyme_list))


def get_stmt_sites(stmts):
    stmt_sites = set([])
    for stmt in stmts:
        site_info = '%s%s' % (stmt.residue, stmt.position)
        stmt_sites.add((stmt.sub.name, site_info))
    return stmt_sites


def coverage(set1, set2):
    coverage = len(set1.intersection(set2))
    return coverage


def load_brca_sites():
    filename = 'sources/Merged_dataset_normalized_subset.csv'
    sites = set([])
    for row in read_unicode_csv(filename, skiprows=1):
        entry_info = row[0]
        site_info = entry_info.split('_')[1]
        up_id = row[-1]
        gene_name = uniprot_client.get_gene_name(up_id)
        sites.add((gene_name, site_info))
    return sites


def get_brca_corr_prior(site_list):
    with open('brca_corr_prot_list.txt', 'rt') as f:
        genes = [line.strip() for line in f]
    prior = {}
    for site in site_list:
        site_key = '%s:%s' % (site[0], site[1])
        prior[site_key] = genes
    return prior


def load_pc_phos():
    site_info = {}
    stmts = []
    for row in read_unicode_csv('sources/PathwayCommons9.All.hgnc.sif',
                                delimiter='\t'):
        relation = row[1]
        if relation != 'controls-phosphorylation-of':
            continue
        enzyme = row[0]
        enz_agent = Agent(enzyme, db_refs={'HGNC': enzyme})
        substrate = row[2]
        sub_agent = Agent(substrate, db_refs={'HGNC': substrate})
        stmt = Phosphorylation(enz_agent, sub_agent)
        stmts.append(stmt)
    return stmts


def save_gene_prior(prior, filename):
    site_prior = {}
    for row in read_unicode_csv('mapped_peptides.txt', delimiter='\t',
                                skiprows=1):
        # Check for uniprot ID
        site_id = row[0]
        gene_name = row[2]
        if gene_name in prior:
            site_prior[site_id] = prior[gene_name]
    with open(filename, 'wt') as f:
        for site_id, gene_list in site_prior.items():
            f.write('%s\t%s\n' % (site_id, ','.join(gene_list)))


def save_default_prior(gene_list, filename):
    with open(filename, 'wt') as f:
        for row in read_unicode_csv('mapped_peptides.txt', delimiter='\t',
                                    skiprows=1):
            # Check for uniprot ID
            site_id = row[0]
            f.write('%s\t%s\n' % (site_id, ','.join(gene_list)))


if __name__ == '__main__':
    #ov_sites = load_ov_sites(use_mapped=True)
    #brca_prior = get_brca_corr_prior(ov_sites)
    #save_prior(brca_prior, 'brca_prot_corr_100_prior.txt')

    #brca_sites = load_brca_sites()
    #print("BRCA phospho-MS data: %d of %d peptides" %
    #      (coverage(ov_sites, brca_sites), len(ov_sites)))

    """
    reg_stmts = get_indra_reg_act_stmts()
    act_stmts = ac.filter_by_type(reg_stmts, Activation)
    inh_stmts = ac.filter_by_type(reg_stmts, Inhibition)
    reg_stmts = act_stmts + inh_stmts
    reg_stmts = [s for s in reg_stmts if s.subj is not None]
    reg_stmts = ac.filter_genes_only(reg_stmts, specific_only=True)
    """

    #indra_stmts = get_indra_phos_stmts()
    """
    indra_stmts = ac.load_statements('sources/indra_phos_stmts.pkl')
    syn_stmts = load_statements_from_synapse(synapse_id='syn10998244')
    pc_stmts = load_pc_phos()
    omni_stmts = get_omnipath_stmts()
    phos_stmts = get_phosphosite_stmts()
    all_stmts = syn_stmts + omni_stmts + phos_stmts + indra_stmts + pc_stmts
    ac.dump_statements(all_stmts, 'sources/all_stmts.pkl')
    """
    all_stmts = ac.load_statements('sources/all_stmts.pkl')
    nsprior = to_nonspec_prior(all_stmts)
    nsprior_filename = 'priors/indra_nkconf2_combined_prot_spec.txt'
    save_gene_prior(nsprior, nsprior_filename)
    syn = synapseclient.login()
    syn_file = synapseclient.File(nsprior_filename, parent='syn11272284')
    syn.store(syn_file)

    all_kinases = [k for kin_list in nsprior.values()
                     for k in kin_list]
    kin_ctr = Counter(all_kinases)
    kin_ctr = sorted([(k, v) for k, v in kin_ctr.items()],
                     key=lambda x: x[1], reverse=True)

    default_prior_list = [t[0] for t in kin_ctr[0:200]]
    default_prior_filename = 'priors/indra_nkconf2_combined_default200.txt'
    save_default_prior(default_prior_list, default_prior_filename)
    syn_file = synapseclient.File(default_prior_filename, parent='syn11272284')
    syn.store(syn_file)

    import sys
    sys.exit()

    """
    with open('sources/stmt_cache.pkl', 'rb') as f:
        syn_stmts, omni_stmts, phos_stmts, indra_stmts = pickle.load(f)
    """

    #db_stmts = syn_stmts + omni_stmts + phos_stmts
    all_stmts = ac.filter_genes_only(all_stmts, specific_only=True)

    all_prior = to_prior(all_stmts)

    ext_prior_100 = add_regulators(reg_stmts, all_prior, max_features=100)
    save_prior(ext_prior_100, 'regulators_prior_100.txt')
    ext_prior_200 = add_regulators(reg_stmts, all_prior, max_features=200)
    save_prior(ext_prior_200, 'regulators_prior_200.txt')


    print("Phosphosite: %d of %d peptides" %
          (coverage(ov_sites, get_stmt_sites(phos_stmts)), len(ov_sites)))
    print("Phosphosite + NetworKIN: %d of %d peptides" %
          (coverage(ov_sites, get_stmt_sites(syn_stmts)), len(ov_sites)))
    print("Omnipath (incl. PSP, Signor, et al.): %d of %d peptides" %
          (coverage(ov_sites, get_stmt_sites(omni_stmts)), len(ov_sites)))
    print("REACH/INDRA: %d of %d peptides" %
          (coverage(ov_sites, get_stmt_sites(indra_stmts)), len(ov_sites)))
    print("Combined prior: %d of %d peptides" %
          (coverage(ov_sites, get_stmt_sites(all_stmts)), len(ov_sites)))
    print("BRCA phospho-MS data: %d of %d peptides" %
          (coverage(ov_sites, brca_sites), len(ov_sites)))

    all_sites = get_stmt_sites(all_stmts).union(brca_sites)
    print("Combined all: %d of %d peptides" %
          (coverage(ov_sites, all_sites), len(ov_sites)))

    db_prior = to_prior(db_stmts)
    # Get activators of kinases
    ext_db_prior = add_regulators(reg_stmts, db_prior)
    ext_prior = add_regulators(reg_stmts, all_prior)

    #indra_prior = to_prior(indra_stmts)
    db_counts = [len(kinases) for kinases in db_prior.values()]
    all_counts = [len(kinases) for kinases in all_prior.values()]
    ext_counts = [len(genes) for genes in ext_prior.values()]
    ext_db_counts = [len(genes) for genes in ext_db_prior.values()]

    plt.ion()
    plt.figure()
    #plt.hist(np.log10(db_counts), bins=20, alpha=0.5)
    #plt.hist(np.log10(all_counts), bins=20, alpha=0.5)
    plt.hist(np.log10(ext_db_counts), bins=20, alpha=0.5)
    plt.hist(np.log10(ext_counts), bins=20, alpha=0.5)

    plt.xlabel('log10(Num annotations)')
    plt.ylabel('Number of peptides')

    # Get list of kinases that bind to the gene


    #save_prior(all_stmts)

    # FOR PLOTTING OVERLAP
    #def get_kin_sub(stmts):
    #    return set([(s.enz.name, s.sub.name, s.position) for s in stmts])
    #indra_sites = get_kin_sub(indra_stmts)
    #phos_sites = get_kin_sub(phos_stmts)
    #nk_sites = get_kin_sub(nk_stmts)


"""
# NOTE: not used anymore because data is being loaded from Synapse
def get_ovarian_nk_stmts():
    stmts = []
    for row_ix, row in enumerate(
                        read_unicode_csv('ovarian_kinase_substrate_table.csv',
                                         skiprows=1)):

        source = row[5]
        sources = set()
        if source != 'NetworKIN':
            sources.add(source)
            continue
        site_info = row[0]
        residue = site_info[0].upper()
        position = site_info[1:]
        enz_hgnc_name = row[2].upper()
        sub_hgnc_name = row[4].upper()
        ev = Evidence(source_api='networkin', source_id='row_%d' % (row_ix+2))
        enz = Agent(enz_hgnc_name, db_refs=get_ids(enz_hgnc_name))
        sub = Agent(sub_hgnc_name, db_refs=get_ids(sub_hgnc_name))
        stmt = Phosphorylation(enz, sub, residue, position, evidence=ev)
        stmts.append(stmt)
    print("Non NK sources: %s" % sources)
    stmts = ac.filter_human_only(stmts)
    stmts = ac.filter_genes_only(stmts)
    return stmts
"""

"""
def get_ovarian_sites():
    ov_sites = []
    for row in read_unicode_csv('ovarian_phosphopeptides.csv',
                                skiprows=1):
        substrate = row[0]
        position = row[2][1:]
        ov_sites.append((substrate, position))
    ov_sites = set(ov_sites)
    return ov_sites
"""

