import sys
import os
os.environ['R_HOME'] = '/nas/longleaf/rhel8/apps/r/4.3.1/lib64/R'
os.environ['R_LIBS_USER'] = '/proj/yunligrp/users/muqing/R-4.3.1'
from scripts.simple import SimpleMethods
from rpy2.robjects import pandas2ri


def main(sample, cluster_num, method):
        ipt_dir = os.path.join("../../../data/DLPFC", sample)
        file_name = '_'.join([sample, "filtered_feature_bc_matrix.h5"])
        tmp_path = '../../../results/tmp/stlearn_tiles'

        methods = SimpleMethods(ipt_dir, file_name, cluster_num, tmp_path)
        clusters = getattr(methods, method)()
        if type(clusters) == 'rpy2.robjects.vectors.ListVector':
            for sample in clusters.names:
                pandas2ri.activate()
                df = pandas2ri.conversion.rpy2py(clusters.rx2(sample))
                df.to_csv(os.path.join("../../../results/tmp/DLPFC/simple", '_'.join([sample, method]) + '.csv'))
        else:
            for sample in clusters:
                df = clusters[sample]
                df.to_csv(os.path.join("../../../results/tmp/DLPFC/simple", '_'.join([sample, method]) + '.csv'))


if __name__ == '__main__':
    sample = sys.argv[1]
    cluster_num = int(sys.argv[2])
    method = sys.argv[3]
    main(sample, cluster_num, method)
    # os.rename(os.path.join(ipt_dir, 'spatial/tissue_positions_list.txt'),
    #           os.path.join(ipt_dir, 'spatial/tissue_positions_list.csv'))
    # os.rename(os.path.join(ipt_dir, 'spatial/tissue_positions_list.csv'),
    #           os.path.join(ipt_dir, 'spatial/tissue_positions_list.txt'))

