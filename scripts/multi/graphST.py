import pandas as pd
import scanpy as sc
from GraphST import GraphST as GST
from GraphST.utils import clustering


def graphST(samples, cnt_paths, xy_paths, genes_paths, cluster_num):
    adata_list = []
    XY = pd.DataFrame()
    for i, sample in enumerate(samples):
        adata = sc.read_mtx(cnt_paths[i])
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        xy['sample_name'] = sample
        XY = pd.concat([XY, xy], axis=0)
        adata.obsm['spatial'] = xy[['row', 'col']].values
        adata.var.index = pd.read_csv(genes_paths[i], index_col=None, header=0)['genes']
        adata.var_names_make_unique()
        adata_list.append(adata)

    adata = sc.concat(adata_list, join='inner')
    adata.obs.reset_index(drop=True, inplace=True)

    model = GST.GraphST(adata)
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
