import stlearn as st
import pandas as pd
import scanpy as sc


def stLearn(samples, cnt_paths, xy_paths, genes_paths, cluster_num, tmp_path):
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

    st.pp.filter_genes(adata, min_cells=1)
    st.pp.normalize_total(adata)
    st.pp.log1p(adata)
    st.em.run_pca(adata, n_comps=15)
    sc.external.pp.harmony_integrate(adata, 'batch')

    XY[['imagecol', 'imagerow']] = XY[['col', 'row']]
    data = st.create_stlearn(count=adata.X, spatial=XY[['imagecol', 'imagerow']], library_id='+'.join(samples))
    data.obsm['X_pca'] = adata.obsm['X_pca_harmony']

    st.pp.tiling(data, tmp_path + "/stlearn_tiles")
    st.pp.extract_feature(data)

    data.obs['array_row'] = data.obs['imagerow']
    data.obs['array_col'] = data.obs['imagecol']
    st.spatial.SME.SME_normalize(data, use_data="raw", weights="physical_distance")
    data.X = data.obsm['raw_SME_normalized']

    st.pp.scale(data)
    st.em.run_pca(data, n_comps=15)

    st.tl.clustering.kmeans(data, n_clusters=cluster_num, use_data="X_pca", key_added="X_harmony_kmeans")

    results = dict()
    clusters = pd.concat([data.obs['X_harmony_kmeans'].reset_index(drop=True), XY[['row', 'col', 'sample_name']].reset_index(drop=True)], axis=1)
    clusters.index = XY.index
    clusters.columns = ['C', 'row', 'col', 'sample_name']
    for sample in samples:
        df = clusters[clusters['sample_name'] == sample]
        del df['sample_name']
        results[sample] = df

    return results

