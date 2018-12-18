import pickle
from collections import Counter

from indra.util import read_unicode_csv
from align_isoforms import load_refseq_up_map
from indra.databases import uniprot_client

import numpy as np
from matplotlib import pyplot as plt
import scipy.stats
import scipy
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import Imputer

import seaborn

from rpy2 import robjects as ro
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()

import synapseclient
syn = synapseclient.login()


cell_lines = [
 'BT20',
 'HCC1806',
 'HS578T',
 'MCF10A',
 'MCF7',
 'MDAMB231',
 'PDX1258',
 'PDX1328',
 'SKBR3',
 'MDAMB134',
 'MDAMB157',
 'MDAMB361',
 'MDAMB436',
 'MDAMB453',
 'MDAMB468',
 'CAL51',
 'CAL851',
 'CAL120',
 'BT549',
 'HCC38',
 'HCC70',
 'HCC1395',
 'HCC1419',
 'HCC1500',
 'HCC1937',
 'HCC1954',
 'PDXHCI002',
 'CAMA1',
 'HCC1143',
 'HCC1428',
 'HME1',
 #'MCF10AREP2',
 'SUM1315',
 'SUM149',
 'SUM159',
 'T47D',
]


def get_mrna_list():
    mrna_set = set()
    for row in read_unicode_csv('sources/gene_types.csv', skiprows=1):
        hgnc_id = row[2]
        gene_type = row[4]
        if gene_type == 'protein_coding':
            mrna_set.add(hgnc_id)
    return list(mrna_set)


def get_cell_line_col_maps(pms_filename, ibaq_filename, rna_filename):
    pms_col_map = {}
    with open(pms_filename, 'rt') as f:
        pms_header = f.readline().strip().split(',')
        for col_ix, col in enumerate(pms_header):
            capcol = col.upper()
            if capcol in cell_lines:
                pms_col_map[capcol] = col_ix

    ibaq_col_map = {}
    with open(ibaq_filename, 'rt') as f:
        ibaq_header = f.readline().strip().split(',')
        for col_ix, col in enumerate(ibaq_header):
            capcol = col.upper()
            if capcol in cell_lines:
                ibaq_col_map[capcol] = col_ix

    rna_col_map = {}
    with open(rna_filename, 'rt') as f:
        rna_header = f.readline().strip().split('\t')
        for col_ix, col in enumerate(rna_header):
            capcol = col.upper()
            if capcol in cell_lines:
                rna_col_map[capcol] = col_ix

    return (pms_col_map, ibaq_col_map, rna_col_map)


def load_data(pms_filename, ibaq_filename, rna_filename):
    pms_col_map, ibaq_col_map, rna_col_map = \
            get_cell_line_col_maps(pms_filename, ibaq_filename, rna_filename)
    # Get the ibaq data first
    prot_labels = []
    prot_rows = []
    for row_ix, row in enumerate(read_unicode_csv(ibaq_filename, skiprows=1)):
        values = []
        gene_name = row[0]
        if gene_name == '0':
            gene_name = None
        uniprot_id = row[1]
        for cell_line in cell_lines:
            col_ix = ibaq_col_map[cell_line]
            if not row[col_ix]:
                val = np.nan
            else:
                val = float(row[col_ix])
            float(val)
            values.append(val)
        values = np.array(values)
        prot_labels.append((gene_name, uniprot_id))
        prot_rows.append(values)

    # Then get the phospho-MS data
    site_labels = []
    site_rows = []
    for row_ix, row in enumerate(read_unicode_csv(pms_filename, skiprows=1)):
        values = []
        site_name = row[0]
        site = site_name.split('_')[1]
        uniprot_id = row[-3]
        uniprot_id_base = row[-1]
        gene_name = uniprot_client.get_gene_name(uniprot_id_base)
        site_position = row[-5]
        site_key = (gene_name, site)
        for cell_line in cell_lines:
            col_ix = pms_col_map[cell_line]
            if not row[col_ix]:
                val = np.nan
            else:
                val = float(row[col_ix])
            float(val)
            values.append(val)
        values = np.array(values)
        site_labels.append((site_key, site_name, uniprot_id, uniprot_id_base,
                            gene_name, site_position))
        site_rows.append(values)

    # Then get the RNA data
    rna_labels = []
    rna_rows = []
    mrna_list = get_mrna_list()
    zero_var = 0
    for row_ix, row in enumerate(read_unicode_csv(rna_filename, skiprows=1,
                                                  delimiter='\t')):
        values = []
        gene_name = row[0]
        if gene_name not in mrna_list:
            continue
        for cell_line in cell_lines:
            col_ix = rna_col_map[cell_line]
            val = float(row[col_ix])
            values.append(val)
        values = np.array(values)
        if np.var(values) == 0:
            zero_var += 1
            continue
        rna_labels.append((gene_name,))
        rna_rows.append(values)
    print("%d RNA entries with zero variance" % zero_var)
    site_arr = np.array(site_rows)
    prot_arr = np.array(prot_rows)
    rna_arr = np.array(rna_rows)

    return {'prot_labels': prot_labels, 'prot_arr': prot_arr,
            'site_labels': site_labels, 'site_arr': site_arr,
            'rna_labels': rna_labels, 'rna_arr': rna_arr}


def get_top_correlations(corrs, site_labels, pred_labels, max_corrs=200):
    site_dict = {}
    for row_ix in range(corrs.shape[0]):
        label = site_labels[row_ix]
        corr_row = corrs[row_ix,:]
        sort_ixs = np.argsort(np.abs(corr_row))
        corr_vec = []
        for ix in sort_ixs[::-1]:
            if len(corr_vec) == max_corrs:
                break
            label = pred_labels[ix]
            if label[0] is not None:
                corr_vec.append((label, corr_row[ix]))
        #corr_vec = [(pred_labels[ix], corr_row[ix])
        #            for ix in sort_ixs[:-max_corrs:-1]]
        site_label = site_labels[row_ix]
        site_dict[site_label] = corr_vec
    return site_dict


def get_default_corrs(corr_dict, num_features=200):
    all_corrs = []
    # Add all the gene names
    for site, corr_list in corr_dict.items():
        all_corrs.extend([t[0][0] for t in corr_list if t[0][0]])

    ctr = Counter(all_corrs)
    ctr = sorted([(k, v) for k, v in ctr.items()], key=lambda x: x[1],
                 reverse=True)
    genes = [t[0] for t in ctr[0:num_features]]
    return genes


def get_site_map(site_labels):
    # For each site in the SC3 phospho data, see if we have a matching site in
    # the BRCA phospho data
    site_map = {}
    brca_site_keys = [t[0] for t in site_labels]
    brca_ix_map = {}
    for ix, brca_site in enumerate(brca_site_keys):
        brca_ix_map[brca_site] = ix
    for row in read_unicode_csv('mapped_peptides.txt', delimiter='\t',
                                skiprows=1):
        site_id = row[0]
        gene_name = row[2]
        orig_site = row[4]
        mapped_site = row[7]
        site_ixs = set()
        for site in ((gene_name, orig_site), (gene_name, mapped_site)):
            brca_ix = brca_ix_map.get(site)
            if brca_ix:
                site_ixs.add(brca_ix)
        # If there's no mapping, don't add to the map
        if not site_ixs:
            continue
        # Otherwise, add
        if site_id in site_map:
            site_map[site_id] |= site_ixs
        else:
            site_map[site_id] = site_ixs
    site_map_list = {}
    for k, v in site_map.items():
        site_map_list[k] = list(v)
    return site_map_list



def build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                prot_default, rna_default, peptide_specific=True,
                num_features=100):
    peptide_file = \
        'sources/retrospective_ova_phospho_sort_common_gene_10057.txt'
    counter = 0
    prior = {}
    for row in read_unicode_csv(peptide_file, delimiter='\t', skiprows=1):
        site_id = row[0]
        gene_sym, rem = site_id.split('.', maxsplit=1)
        rs_id, site_info = rem.split(':')
        if site_id in site_map and peptide_specific:
            brca_site_ix_list = site_map[site_id]
            if len(brca_site_ix_list) > 1:
                print("More than one site for %s" % site_id)
            brca_site_ix = brca_site_ix_list[0] # FIXME
            brca_site = site_labels[brca_site_ix]
            prot_prior = [t[0][0] for t in prot_corr_dict[brca_site]]
            rna_prior = [t[0][0] for t in rna_corr_dict[brca_site]]
            prior[site_id] = (prot_prior[0:num_features],
                              rna_prior[0:num_features])
            counter += 1
        else:
            prior[site_id] = (prot_default[0:num_features],
                              rna_default[0:num_features])
    print("%d peptide-specific priors used" % counter)
    return prior


def save_prior(filename, prior, prot_only=False):
    with open(filename, 'wt') as f:
        for site_id, (prot_prior, rna_prior) in prior.items():
            prot_prior_str = ','.join(prot_prior)
            rna_prior_str = ','.join(rna_prior)
            if prot_only:
                line = '%s\t%s\n' % (site_id, prot_prior_str)
            else:
                line = '%s\t%s\t%s\n' % (site_id, rna_prior_str, prot_prior_str)
            f.write(line)


def vip(x, y, model):
    t = model.x_scores_
    w = model.x_weights_
    q = model.y_loadings_

    m, p = x.shape
    _, h = t.shape

    vips = np.zeros((p,))

    s = np.diag(t.T @ t @ q.T @ q).reshape(h, -1)
    total_s = np.sum(s)

    for i in range(p):
        weight = np.array([ (w[i,j] / scipy.linalg.norm(w[:,j]))**2
                            for j in range(h) ])
        #weight = np.array([ (w[i,j] / np.linalg.norm(w[:,j]))**2
        #                    for j in range(h) ])
        vips[i] = np.sqrt(p*(s.T @ weight)/total_s)
    return vips


def prior_from_scores(scores, labels, num_features=100):
    sorted_scores = np.argsort(scores)
    prior = set()
    ind = []
    for pred_ix in sorted_scores[::-1]:
        if len(prior) == num_features:
            break
        pred_label = labels[pred_ix]
        pred_gene = pred_label[0]
        if not pred_gene:
            continue
        else:
            prior.add(pred_gene)
    return list(prior)


def save_and_upload(filename, prior):
    save_prior(filename, prior)
    syn_file = synapseclient.File(filename, parent='syn11272284')
    syn.store(syn_file)


if __name__ == '__main__':
    pms_filename = 'sources/Merged_dataset_normalized_subset.csv'
    ibaq_filename = 'sources/ibaq_normalized.csv'
    rna_filename = 'sources/RNAseq-rpkm.tsv'

    print("Loading data")
    data = load_data(pms_filename, ibaq_filename, rna_filename)
    site_labels = data['site_labels']
    site_map = get_site_map(site_labels)

    """
    print("Clustering")
    X = data['site_arr']
    Z = linkage(X, method='ward')
    dendrogram(Z, leaf_rotation=90., leaf_font_size=8.)
    plt.show()
    #c, coph_dists = cophenet(Z, pdist(X, metric='correlation'))
    #print(c)
    res = seaborn.clustermap(data=X, z_score=0)
    plt.show()
    """

    print("Imputing values")
    imp = Imputer()
    sa = imp.fit_transform(data['site_arr'].T)
    pa = imp.fit_transform(data['prot_arr'].T)
    ra = imp.fit_transform(data['rna_arr'].T)

    # Calculate the correlations
    recalculate_corrs = False
    if recalculate_corrs:
        print("Calculating protein correlation coefficients")
        prot_corrs = np.array(ro.r.cor(sa, pa, method='spearman'))
        np.save('brca_spearman_prot', prot_corrs)

        print("Calculating RNA correlation coefficients")
        rna_corrs = np.array(ro.r.cor(sa, ra, method='spearman'))
        np.save('brca_spearman_rna', rna_corrs)
    else:
        print("Loading correlations")
        prot_corrs = np.load('brca_spearman_prot.npy')
        rna_corrs = np.load('brca_spearman_rna.npy')

    do_cca = True
    if do_cca:
        print("Running CCA on protein data")
        prot_cca = CCA(n_components=33, scale=False)
        prot_cca.fit(pa, sa)

        print("Running CCA on RNA data")
        rna_cca = CCA(n_components=33, scale=False)
        rna_cca.fit(ra, sa)

        prot_vips = vip(pa, sa, prot_cca)
        prot_vip_prior = prior_from_scores(prot_vips, data['prot_labels'])
        prot_wt_prior = prior_from_scores(np.mean(prot_cca.x_weights_, axis=1),
                                          data['prot_labels'])
        prot_ld_prior = prior_from_scores(np.mean(prot_cca.x_loadings_, axis=1),
                                          data['prot_labels'])

        rna_vips = vip(ra, sa, rna_cca)
        rna_vip_prior = prior_from_scores(rna_vips, data['rna_labels'])
        rna_wt_prior = prior_from_scores(np.mean(rna_cca.x_weights_, axis=1),
                                         data['rna_labels'])
        rna_ld_prior = prior_from_scores(np.mean(rna_cca.x_loadings_, axis=1),
                                         data['rna_labels'])

        # Build the CCA-VIP prior
        prior = build_prior(site_map, site_labels, None,
                            None, prot_vip_prior, rna_vip_prior,
                            peptide_specific=False, num_features=100)
        filename = 'priors/brca_prot_rna_cca_vips_corr100.txt'
        save_and_upload(filename, prior)

        # Build the CCA-Weights prior
        prior = build_prior(site_map, site_labels, None,
                            None, prot_wt_prior, rna_wt_prior,
                            peptide_specific=False, num_features=100)
        filename = 'priors/brca_prot_rna_cca_avgwts_corr100.txt'
        save_and_upload(filename, prior)

        # Build the CCA-Loadings prior
        prior = build_prior(site_map, site_labels, None,
                            None, prot_ld_prior, rna_ld_prior,
                            peptide_specific=False, num_features=100)
        filename = 'priors/brca_prot_rna_cca_avgloadings_corr100.txt'
        save_and_upload(filename, prior)

    print("Getting top correlations")
    # Get the top correlations for each site
    prot_corr_dict = get_top_correlations(prot_corrs,
                                     data['site_labels'], data['prot_labels'],
                                     max_corrs=400)

    rna_corr_dict = get_top_correlations(rna_corrs,
                                     data['site_labels'], data['rna_labels'],
                                     max_corrs=400)

    # Get the default priors
    prot_default = get_default_corrs(prot_corr_dict, num_features=400)
    rna_default = get_default_corrs(rna_corr_dict, num_features=400)

    # Build specific RNA+protein prior with 50 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=True,
                        num_features=50)
    filename = 'priors/brca_prot_rna_specific_corr50.txt'
    save_and_upload(filename, prior)

    # Build the peptide-specific prior
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=True,
                        num_features=100)
    filename = 'priors/brca_prot_rna_specific_corr100.txt'
    save_and_upload(filename, prior)

    # Build the peptide-specific prior with 200 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=True,
                        num_features=200)
    filename = 'priors/brca_prot_rna_specific_corr200.txt'
    save_and_upload(filename, prior)

    # Build the peptide-specific prior with 400 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=True,
                        num_features=400)
    filename = 'priors/brca_prot_rna_specific_corr400.txt'
    save_and_upload(filename, prior)

    # Build the default RNA+protein prior with 100 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=False,
                        num_features=100)
    filename = 'priors/brca_prot_rna_default_corr100.txt'
    save_and_upload(filename, prior)

    # Build protein-only default prior with 200 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=False,
                        num_features=200)
    filename = 'priors/brca_prot_rna_default_corr200.txt'
    save_and_upload(filename, prior)

    # Build protein-only default prior with 200 features
    prior = build_prior(site_map, site_labels, prot_corr_dict, rna_corr_dict,
                        prot_default, rna_default, peptide_specific=False,
                        num_features=400)
    filename = 'priors/brca_prot_rna_default_corr400.txt'
    save_and_upload(filename, prior)

