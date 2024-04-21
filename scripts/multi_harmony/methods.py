import os
import rpy2.robjects as robjects

from .graphST import graphST as run_graphST_multi_harmony
from .spaGCN import spaGCN as run_spaGCN_multi_harmony
from .stLearn import stLearn as run_stLearn_multi_harmony
from .STAGATE import STAGATE as run_STAGATE_multi_harmony

scripts_dir = os.path.dirname(os.path.realpath(__file__))
robjects.r.source(os.path.join(scripts_dir, 'BASS.R'))
robjects.r.source(os.path.join(scripts_dir, 'BayesSpace.R'))
robjects.r.source(os.path.join(scripts_dir, 'MAPLE.R'))
robjects.r.source(os.path.join(scripts_dir, 'Seurat.R'))
robjects.r.source(os.path.join(scripts_dir, 'Stardust.R'))


class MultiHarmonyMethods:
    def __init__(self, samples, cnt_paths, xy_paths, genes_paths, img_paths, tmp_path, cluster_num):
        self.__dict__.update(locals())
        del self.__dict__['self']

    def Bass(self):
        return robjects.r.run_Bass_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def BayesSpace(self):
        return robjects.r.run_BayesSpace_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def GraphST(self):
        return run_graphST_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def STAGATE(self):
        return run_STAGATE_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def MAPLE(self):
        return robjects.r.run_MAPLE_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def SpaGCN(self):
        return run_spaGCN_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.img_paths, self.cluster_num)

    def Seurat(self):
        return robjects.r.run_Seurat_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def Stardust(self):
        return robjects.r.run_Sardust_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num)

    def stlearn(self):
        return run_stLearn_multi_harmony(self.samples, self.cnt_paths, self.xy_paths, self.genes_paths, self.cluster_num, self.tmp_path)


if __name__ == '__main__':
    pass
