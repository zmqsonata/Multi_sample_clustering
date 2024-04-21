import pandas as pd
import scanpy as sc
import STAGATE_pyG


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
    adata.uns['Spatial_Net'] = pd.concat([data.uns['Spatial_Net'] for data in adata_list])

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata = STAGATE_pyG.train_STAGATE(adata, n_epochs=500)
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
