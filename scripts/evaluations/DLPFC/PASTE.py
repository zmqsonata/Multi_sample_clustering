import os
import scanpy as sc
import paste as pst
import pandas as pd
from scipy.io import mmwrite


def paste_align(samples, cnt_paths, xy_paths, genes_paths, opt_dir):
    slices = []
    ori_slices = []
    xy_lists = []
    for i, sample in enumerate(samples):
        adata = sc.read_mtx(cnt_paths[i])
        xy = pd.read_csv(xy_paths[i], index_col=0, header=0)
        xy_lists.append(xy)
        adata.obsm['spatial'] = xy[['row', 'col']].values
        adata.var.index = pd.read_csv(genes_paths[i], index_col=None, header=0)['genes']
        adata.var_names_make_unique()
        adata.obs.index = xy.index
        ori_slices.append(adata)
        sc.pp.filter_genes(adata, min_counts=15)
        sc.pp.filter_cells(adata, min_counts=100)
        slices.append(adata)

    pst.filter_for_common_genes(slices)

    initial_slice = slices[0].copy()
    lmbda = len(slices) * [1 / len(slices)]

    b = [pst.match_spots_using_spatial_heuristic(slices[0].X.toarray(), slice.X.toarray()) for slice in slices]
    center_slice, pis = pst.center_align(initial_slice, slices, lmbda, random_seed=5, pis_init=b)
    _, new_slices = pst.stack_slices_center(center_slice, slices, pis)

    for i, sample in enumerate(samples):
        new_xy = xy_lists[i]
        new_xy = new_xy.loc[new_slices[i].obs.index]
        new_xy[['row', 'col']] = new_slices[i].obsm['spatial']
        new_cnts = ori_slices[i][new_xy.index, :]
        mmwrite(os.path.join(opt_dir, sample + "_cnt_aligned.mtx"), new_cnts.X)
        new_xy.to_csv(os.path.join(opt_dir, sample + "_xy_aligned.csv"))


if __name__ == "__main__":
    samples_nums = (("151507", "151508", "151509", "151510"),
                    ("151669", "151670", "151671", "151672"), ("151673", "151674", "151675", "151676"))
    dir = "../../../results/tmp"
    for samples in samples_nums:
        cnt_paths = []
        xy_paths = []
        genes_paths = []
        for sample in samples:
            cnt_paths.append(os.path.join(dir, sample + "_cnt.mtx"))
            xy_paths.append(os.path.join(dir, sample + "_xy.csv"))
            genes_paths.append(os.path.join(dir, sample + "_genes.csv"))
        paste_align(samples, cnt_paths, xy_paths, genes_paths, dir)
