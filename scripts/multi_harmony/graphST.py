import pandas as pd
import scanpy as sc
from GraphST import GraphST as GST
from GraphST.utils import clustering
import anndata


def graphST(samples, cnt_paths, xy_paths, genes_paths, cluster_num):
    adata_list = []
    XY = pd.DataFrame()
    for i, sample in enumerate(samples):
        adata = sc.read_mtx(cnt_paths[i])
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        xy['sample_name'] = sample
        XY = pd.concat([XY, xy], axis=0)
        adata.var.index = pd.read_csv(genes_paths[i], index_col=None, header=0)['genes']
        adata.var_names_make_unique()
        adata_list.append(adata)

    adata = sc.concat(adata_list, join='inner')
    adata.obs.reset_index(drop=True, inplace=True)
    adata.obs['batch'] = XY['sample_name'].reset_index(drop=True)

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata, n_comps=50)
    sc.external.pp.harmony_integrate(adata, 'batch')

    adata_new = anndata.AnnData(X=adata.obsm['X_pca_harmony'])
    adata_new.var['highly_variable'] = True
    adata_new.obsm['spatial'] = XY[['row', 'col']].values
    model = GST.GraphST(adata_new)
    adata = model.train()

    tool = 'mclust'
    if tool == 'mclust':
        clustering(adata, cluster_num, method=tool)
    elif tool in ['leiden', 'louvain']:
        clustering(adata, cluster_num, method=tool, start=0.1, end=2.0, increment=0.01)

    results = dict()
    clusters = pd.concat([adata.obs['mclust'].reset_index(drop=True), XY[['row', 'col', 'sample_name']].reset_index(drop=True)], axis=1)
    clusters.index = XY.index
    clusters.columns = ['C', 'row', 'col', 'sample_name']
    for sample in samples:
        df = clusters[clusters['sample_name'] == sample]
        del df['sample_name']
        results[sample] = df

    return results

