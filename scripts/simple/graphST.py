import torch
import scanpy as sc
from GraphST.utils import clustering
from GraphST import GraphST


def run_graphST_Simple(dir, file_name, cluster_num, use_gpu=True):
    use_gpu = True if torch.cuda.is_available() and use_gpu else False
    device = torch.device('cuda:0' if use_gpu else 'cpu')
    adata = sc.read_visium(dir, count_file=file_name, load_images=True)
    adata.var_names_make_unique()

    model = GraphST.GraphST(adata, device=device)

    adata = model.train()

    radius = 50

    tool = 'mclust'

    # clustering
    if tool == 'mclust':
       clustering(adata, cluster_num, radius=radius, method=tool, refinement=True)
    elif tool in ['leiden', 'louvain']:
       clustering(adata, cluster_num, radius=radius, method=tool, start=0.1, end=2.0, increment=0.01, refinement=False)

    results = adata.obs[['mclust', 'array_row', 'array_col']]
    results.columns = ['C', 'row', 'col']
    return results
