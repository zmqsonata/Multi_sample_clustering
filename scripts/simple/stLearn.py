import stlearn as st
from pathlib import Path


def run_stLearn_Simple(dir, file_name, cluster_num, tmp_path):

    TILE_PATH = Path(tmp_path)/"stlearn_tiles"
    TILE_PATH.mkdir(parents=True, exist_ok=True)
    data = st.Read10X(dir, count_file=file_name)

    st.pp.filter_genes(data,min_cells=1)
    st.pp.normalize_total(data)
    st.pp.log1p(data)

    st.em.run_pca(data,n_comps=15)

    st.pp.tiling(data, TILE_PATH)

    st.pp.extract_feature(data)

    # stSME
    st.spatial.SME.SME_normalize(data, use_data="raw", weights="physical_distance")
    data_ = data.copy()
    data_.X = data_.obsm['raw_SME_normalized']

    st.pp.scale(data_)
    st.em.run_pca(data_,n_comps=15)

    st.tl.clustering.kmeans(data_, n_clusters=cluster_num, use_data="X_pca", key_added="X_pca_kmeans")

    results = data_.obs[['X_pca_kmeans', 'array_row', 'array_col']]
    results.columns = ['C', 'row', 'col']
    return results
