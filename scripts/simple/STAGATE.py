import warnings
warnings.filterwarnings("ignore")
import scanpy as sc
import STAGATE_pyG

def run_STAGATE_Simple(dir, file_name, cluster_num):
    adata = sc.read_visium(path=dir, count_file=file_name)
    adata.var_names_make_unique()

    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=3000)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    STAGATE_pyG.Cal_Spatial_Net(adata, rad_cutoff=150)
    adata = STAGATE_pyG.train_STAGATE(adata)

    sc.pp.neighbors(adata, use_rep='STAGATE')
    sc.tl.umap(adata)

    adata = STAGATE_pyG.mclust_R(adata, used_obsm='STAGATE', num_cluster=cluster_num)
    results = adata.obs[['mclust', 'array_row', 'array_col']]
    results.columns = ['C', 'row', 'col']
    return results
