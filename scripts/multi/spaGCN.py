import pandas as pd
import scanpy as sc
import random
import torch
import SpaGCN as spg
from PIL import Image
import numpy as np
from SpaGCN.calculate_adj import pairwise_distance


def calculate_adj_matrix(x, y, x_pixel=None, y_pixel=None, images=None, beta=49, alpha=1):
    print("Calculating adj matrix using histology image...")
    beta_half =round(beta/2)
    g=[]
    for i in range(len(x_pixel)):
        for j in range(len(x_pixel[i])):
            max_x=images[i].shape[0]
            max_y=images[i].shape[1]
            nbs=images[i][max(0,x_pixel[i][j]-beta_half):min(max_x,x_pixel[i][j]+beta_half+1),
                max(0,y_pixel[i][j]-beta_half):min(max_y,y_pixel[i][j]+beta_half+1)]
            g.append(np.mean(np.mean(nbs,axis=0),axis=0))
    c0, c1, c2=[], [], []
    for i in g:
        c0.append(i[0])
        c1.append(i[1])
        c2.append(i[2])
    c0=np.array(c0)
    c1=np.array(c1)
    c2=np.array(c2)
    print("Var of c0,c1,c2 = ", np.var(c0),np.var(c1),np.var(c2))
    c3=(c0*np.var(c0)+c1*np.var(c1)+c2*np.var(c2))/(np.var(c0)+np.var(c1)+np.var(c2))
    c4=(c3-np.mean(c3))/np.std(c3)
    z_scale=np.max([np.std(x), np.std(y)])*alpha
    z=c4*z_scale
    z=z.tolist()
    print("Var of x,y,z = ", np.var(x),np.var(y),np.var(z))
    X=np.array([x, y, z]).T.astype(np.float32)
    return pairwise_distance(X)


def spaGCN(samples, cnt_paths, xy_paths, genes_paths, img_paths, cluster_num):
    adata_list = []
    img_list = []
    xy_list = []
    for i, sample in enumerate(samples):
        adata = sc.read_mtx(cnt_paths[i])
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        xy_list.append(xy)
        xy['sample_name'] = sample
        adata.var.index = pd.read_csv(genes_paths[i], index_col=None, header=0)['genes']
        adata.var_names_make_unique()
        adata_list.append(adata)
        img_list.append(np.array(Image.open(img_paths[i])))

    XY = pd.concat(xy_list, axis=0)
    adata = sc.concat(adata_list, join='inner')
    adata.obs.reset_index(drop=True, inplace=True)

    s = 1
    b = 49
    x_pixels = [xy['x_pixel'].tolist() for xy in xy_list]
    y_pixels = [xy['y_pixel'].tolist() for xy in xy_list]
    adj = calculate_adj_matrix(x=XY['row'], y=XY['col'], x_pixel=x_pixels, y_pixel=y_pixels, images=img_list, beta=b, alpha=s)

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

    clf.train(adata, adj, init_spa=True, init="louvain", res=res, tol=5e-3, lr=0.05, max_epochs=20000)
    y_pred, prob = clf.predict()
    adata.obs["pred"] = y_pred
    adata.obs["pred"] = adata.obs["pred"].astype('category')

    # shape="hexagon" for Visium data, "square" for ST data.
    adj_2d = spg.calculate_adj_matrix(x=XY['row'], y=XY['col'], histology=False)
    refined_pred = spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d,
                              shape="hexagon")
    adata.obs["pred"] = refined_pred

    results = dict()
    clusters = pd.concat([adata.obs["pred"].reset_index(drop=True), XY[['row', 'col', 'sample_name']].reset_index(drop=True)], axis=1)
    clusters.index = XY.index
    clusters.columns = ['C', 'row', 'col', 'sample_name']
    for sample in samples:
        df = clusters[clusters['sample_name'] == sample]
        del df['sample_name']
        results[sample] = df

    return results
