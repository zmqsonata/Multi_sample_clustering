from scipy.io import mmread
import stlearn as st
import pandas as pd


def stLearn(samples, cnt_paths, xy_paths, genes_paths, cluster_num, tmp_path):
    cnt_list = []
    genes_list = []
    XY = pd.DataFrame()
    for i, sample in enumerate(samples):
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        xy['sample_name'] = sample
        XY = pd.concat([XY, xy], axis=0)
        genes_list.append(pd.read_csv(genes_paths[i], index_col=None, header=0)['genes'])

    common_genes = list(set.intersection(*[set(genes) for genes in genes_list]))
    for i in range(len(samples)):
        mtx = pd.DataFrame.sparse.from_spmatrix(mmread(cnt_paths[i]))
        mtx.columns = genes_list[i]
        cnt_list.append(mtx[common_genes])

    cnt = pd.concat(cnt_list)
    cnt.index = XY.index
    XY[['imagecol', 'imagerow']] = XY[['col', 'row']]
    data = st.create_stlearn(count=cnt, spatial=XY[['imagecol', 'imagerow']], library_id='+'.join(samples))

    st.pp.filter_genes(data,min_cells=1)
    st.pp.normalize_total(data)
    st.pp.log1p(data)
    st.em.run_pca(data, n_comps=15)
    st.pp.tiling(data, tmp_path + "/stlearn_tiles")
    st.pp.extract_feature(data)

    data.obs['array_row'] = data.obs['imagerow']
    data.obs['array_col'] = data.obs['imagecol']
    st.spatial.SME.SME_normalize(data, use_data="raw", weights="physical_distance")
    data.X = data.obsm['raw_SME_normalized']

    st.pp.scale(data)
    st.em.run_pca(data, n_comps=15)
    st.tl.clustering.kmeans(data, n_clusters=cluster_num, use_data="X_pca", key_added="X_pca_kmeans")

    results = dict()
    clusters = pd.concat([data.obs['X_pca_kmeans'].reset_index(drop=True), XY[['row', 'col', 'sample_name']].reset_index(drop=True)], axis=1)
    clusters.index = XY.index
    clusters.columns = ['C', 'row', 'col', 'sample_name']
    for sample in samples:
        df = clusters[clusters['sample_name'] == sample]
        del df['sample_name']
        results[sample] = df

    return results
