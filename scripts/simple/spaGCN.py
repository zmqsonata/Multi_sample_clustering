import numpy as np
import random
import torch
import warnings
import SpaGCN as spg
import scanpy as sc
warnings.filterwarnings("ignore")


def run_spaGCN_Simple(dir, file_name, cluster_num):
    adata = sc.read_visium(dir, count_file=file_name, load_images=True)

    img = next(iter(adata.uns['spatial'].values()))
    pixels = img['images']['hires']
    scale = img['scalefactors']['tissue_hires_scalef']

    adata.obs[["x_pixel", "y_pixel"]] = adata.obsm['spatial']*scale
    adata.obs[["x_pixel", "y_pixel"]] = adata.obs[["x_pixel", "y_pixel"]].astype('int')
    x_pixel = adata.obs["x_pixel"].tolist()
    y_pixel = adata.obs["y_pixel"].tolist()
    adata = adata[adata.obs["in_tissue"] == 1]
    adata.var_names = [i.upper() for i in list(adata.var_names)]
    adata.var["genename"] = adata.var.index.astype("str")
    adata.var_names_make_unique()

    s = 1
    b = 49
    adj = spg.calculate_adj_matrix(x=x_pixel, y=y_pixel, x_pixel=x_pixel, y_pixel=y_pixel, image=pixels, beta=b, alpha=s,
                                   histology=True)

    spg.prefilter_genes(adata, min_cells=3)
    spg.prefilter_specialgenes(adata)
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)

    p = 0.5
    l = spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

    r_seed = t_seed = n_seed = 100

    res = spg.search_res(adata, adj, l, cluster_num, start=0.7, step=0.1, tol=5e-3, lr=0.05, max_epochs=20,
                         r_seed=r_seed,
                         t_seed=t_seed, n_seed=n_seed)

    clf = spg.SpaGCN()
    clf.set_l(l)

    random.seed(r_seed)
    torch.manual_seed(t_seed)
    np.random.seed(n_seed)

    clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')
    # Do cluster refinement(optional)
    # shape="hexagon" for Visium data, "square" for ST data.
    adj_2d = spg.calculate_adj_matrix(x=adata.obs['array_row'], y=adata.obs['array_col'], histology=False)
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d,
                              shape="hexagon")
    adata.obs["refined_pred"] = refined_pred
    results = adata.obs[['refined_pred', 'array_row', 'array_col']]
    results.columns = ['C', 'row', 'col']
    return results
