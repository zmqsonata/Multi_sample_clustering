import pandas as pd
import scanpy as sc
import STAGATE_pyG
import anndata

import os
os.environ['R_LIBS_USER'] = '/proj/yunligrp/users/muqing/R-4.3.1'
os.environ['R_HOME'] = '/nas/longleaf/rhel8/apps/r/4.3.1/lib64/R'


def STAGATE(samples, cnt_paths, xy_paths, genes_paths, cluster_num):
    adata_list = []
    XY = pd.DataFrame()
    for i, sample in enumerate(samples):
        adata = sc.read_mtx(cnt_paths[i])
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        adata.obsm['spatial'] = xy[['row', 'col']].values
        xy['sample_name'] = sample
        XY = pd.concat([XY, xy], axis=0)
        adata.var.index = pd.read_csv(genes_paths[i], index_col=None, header=0)['genes']
        adata.var_names_make_unique()
        sc.pp.calculate_qc_metrics(adata, inplace=True)
        STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=11)
        adata_list.append(adata)

    adata = sc.concat(adata_list, join='inner')
    adata.obs.reset_index(drop=True, inplace=True)
    adata.obs['batch'] = XY['sample_name'].reset_index(drop=True)
    adata.uns['Spatial_Net'] = pd.concat([data.uns['Spatial_Net'] for data in adata_list])

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=50)
    sc.external.pp.harmony_integrate(adata, 'batch')

    adata_new = anndata.AnnData(X=adata.obsm['X_pca_harmony'])
    adata_new.uns['Spatial_Net'] = adata.uns['Spatial_Net']

    adata = STAGATE_pyG.train_STAGATE(adata_new, n_epochs=500)
    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)
    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=cluster_num)

    results = dict()
    clusters = pd.concat(
        [adata.obs["mclust"].reset_index(drop=True), XY[['row', 'col', 'sample_name']].reset_index(drop=True)], axis=1)
    clusters.index = XY.index
    clusters.columns = ['C', 'row', 'col', 'sample_name']
    for sample in samples:
        df = clusters[clusters['sample_name'] == sample]
        del df['sample_name']
        results[sample] = df

    return results


if __name__ == '__main__':
    cnt_paths = ("../../results/tmp/151507_cnt.mtx", "../../results/tmp/151508_cnt.mtx")
    xy_paths = ("../../results/tmp/151507_xy.csv", "../../results/tmp/151508_xy.csv")
    genes_paths = ("../../results/tmp/151507_genes.csv", "../../results/tmp/151508_genes.csv")
    samples = ('151507', '151508')
    img_paths = ("../../results/tmp/151507_img.png", "../../results/tmp/151508_img.png")
    cluster_num = 5
    STAGATE(samples, cnt_paths, xy_paths, genes_paths, cluster_num)
