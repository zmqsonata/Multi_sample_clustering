import sys
import os
os.environ['R_HOME'] = '/nas/longleaf/rhel8/apps/r/4.3.1/lib64/R'
os.environ['R_LIBS_USER'] = '/proj/yunligrp/users/muqing/R-4.3.1'
from scripts.multi import MultiMethods
from rpy2.robjects import pandas2ri


def main(sample, cluster_num, method):
    img_paths = []
    tmp_path = '../../../results/tmp/stlearn_tiles'
    cnt_paths = ["../../../results/tmp/simulation2/" + sample + "_cnt.mtx"]
    xy_paths = ["../../../results/tmp/simulation2/" + sample + "_xy.csv"]
    genes_paths = ["../../../results/tmp/simulation2/" + sample + "_genes.csv"]

    methods = MultiMethods(sample, cnt_paths, xy_paths, genes_paths, img_paths, tmp_path, cluster_num)
    clusters = getattr(methods, method)()
    if type(clusters) == 'rpy2.robjects.vectors.ListVector':
        for sample in clusters.names:
            pandas2ri.activate()
            df = pandas2ri.conversion.rpy2py(clusters.rx2(sample))
            df.to_csv(os.path.join("../../../results/tmp/simulation2/simple", '_'.join([sample, method]) + '.csv'))
    else:
        for sample in clusters:
            df = clusters[sample]
            df.to_csv(os.path.join("../../../results/tmp/simulation2/simple", '_'.join([sample, method]) + '.csv'))


if __name__ == '__main__':
    sample = sys.argv[1]
    cluster_num = int(sys.argv[2])
    method = sys.argv[3]
    main(sample, cluster_num, method)
