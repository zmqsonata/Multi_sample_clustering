import os
import rpy2.robjects as robjects

from .graphST import run_graphST_Simple
from .STAGATE import run_STAGATE_Simple
from .spaGCN import run_spaGCN_Simple
from .stLearn import run_stLearn_Simple

scripts_dir = os.path.dirname(os.path.realpath(__file__))
robjects.r.source(os.path.join(scripts_dir, 'BASS.R'))
robjects.r.source(os.path.join(scripts_dir, 'BayesSpace.R'))
robjects.r.source(os.path.join(scripts_dir, 'ISC_MEB.R'))
robjects.r.source(os.path.join(scripts_dir, 'MAPLE.R'))
robjects.r.source(os.path.join(scripts_dir, 'Seurat.R'))
robjects.r.source(os.path.join(scripts_dir, 'Stardust.R'))


class SimpleMethods:
    def __init__(self, ipt_dir, file_name, cluster_num, tmp_path):
        self.ipt_dir = ipt_dir
        self.file_name = file_name
        self.cluster_num = cluster_num
        self.tmp_path = tmp_path

    def Bass(self):
        return robjects.r.run_Bass_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def BayesSpace(self):
        return robjects.r.run_BayesSpace_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def GraphST(self):
        return run_graphST_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def STAGATE(self):
        return run_STAGATE_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def ISC_MEB(self):
        return robjects.r.run_ISC_MEB_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def MAPLE(self):
        return robjects.r.run_MAPLE_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def SpaGCN(self):
        return run_spaGCN_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def Seurat(self):
        return robjects.r.run_Seurat_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def Stardust(self):
        return robjects.r.run_Sardust_Simple(self.ipt_dir, self.file_name, self.cluster_num)

    def stlearn(self):
        return run_stLearn_Simple(self.ipt_dir, self.file_name, self.cluster_num, self.tmp_path)


if __name__ == '__main__':
    pass