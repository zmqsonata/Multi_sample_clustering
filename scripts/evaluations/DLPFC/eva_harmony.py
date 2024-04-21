import sys
import os
os.environ['R_HOME'] = '/nas/longleaf/rhel8/apps/r/4.3.1/lib64/R'
os.environ['R_LIBS_USER'] = '/proj/yunligrp/users/muqing/R-4.3.1'
from scripts.multi_harmony import MultiHarmonyMethods
from rpy2.robjects import pandas2ri


def main(samples, cluster_num, method):
    cnt_paths = []
    xy_paths = []
    genes_paths = []
    img_paths = []
    tmp_path = '../../../results/tmp/stlearn_tiles'
    for sample in samples:
        cnt_paths.append("../../../results/tmp/DLPFC/" + sample + "_cnt.mtx")
        xy_paths.append("../../../results/tmp/DLPFC/" + sample + "_xy.csv")
        genes_paths.append("../../../results/tmp/DLPFC/" + sample + "_genes.csv")
        img_paths.append("../../../results/tmp/DLPFC/" + sample + "_img.png")

    methods = MultiHarmonyMethods(samples, cnt_paths, xy_paths, genes_paths, img_paths, tmp_path, cluster_num)
    clusters = getattr(methods, method)()
    if type(clusters) == 'rpy2.robjects.vectors.ListVector':
        for sample in clusters.names:
            pandas2ri.activate()
            df = pandas2ri.conversion.rpy2py(clusters.rx2(sample))
            df.to_csv(os.path.join("../../../results/tmp/DLPFC/harmony", '_'.join([sample, method]) + '.csv'))
    else:
        for sample in clusters:
            df = clusters[sample]
            df.to_csv(os.path.join("../../../results/tmp/DLPFC/harmony", '_'.join([sample, method]) + '.csv'))


if __name__ == '__main__':
    sample = sys.argv[1].split(',')
    cluster_num = int(sys.argv[2])
    method = sys.argv[3]
    main(sample, cluster_num, method)
