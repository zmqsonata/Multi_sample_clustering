# Multi-Slide Clustering

We provide a comprehensive summary and comparison of clustering methods for multi-slide ST data. We outline ten state-of-art ST clustering methods, including seven single-slide ST clustering methods and three multi-slide ST clustering methods, emphasizing crucial aspects such as the statistical or machine learning methodology employed, and the unique features specific to each method. We assess the performance of these methods using two simulation datasets and three real ST datasets. Additionally, we evaluate the effectiveness of data pre-processing methods, including Harmony for gene expression harmonization and PASTE for spatial coordinate alignment, for multi-slide data.



## Methods Compared

|   Method   | Language |                         URL                          | Reference |
| :--------: | :------: | :--------------------------------------------------: | :-------: |
|    BASS    |    R     |          https://github.com/zhengli09/BASS           |    [1]    |
| BayesSpace |    R     |      https://github.com/edward130603/BayesSpace      |    [2]    |
|  iSC_MEB   |    R     |       https://github.com/XiaoZhangryy/iSC.MEB        |    [3]    |
|   MAPLE    |    R     |        https://github.com/carter-allen/maple         |    [4]    |
|   Seurat   |    R     |         https://github.com/satijalab/seurat          |    [5]    |
|  Stardust  |    R     |        https://github.com/InfOmics/stardust/         |    [6]    |
|  GraphST   |  python  |      https://github.com/JinmiaoChenLab/GraphST       |    [7]    |
|   SpaGCN   |  python  |        https://github.com/jianhuupenn/SpaGCN         |    [8]    |
|  STAGATE   |  python  |       https://github.com/QIFEIDKN/STAGATE_pyG        |    [9]    |
|  stLearn   |  python  | https://github.com/BiomedicalMachineLearning/stLearn |    [10]   |

## Repo Structure

- `data`: simulation data
- `results`: results for tests on all datasets
- `scripts`: python and R scripts for running the methods and python scripts for evaluations on DLPFC and simulations datasets
- `visualizations`: Jupyter notebook files for visualize all results on DLPFC and simulations datasets

## Reference
- [1] Zheng Li and Xiang Zhou. “BASS: multi-scale and multi-sample analysis enables accurate cell type clustering and spatial domain detection in spatial transcriptomic studies”. In: Genome biology 23.1 (2022), p. 168.
- [2] Edward Zhao et al. “Spatial transcriptomics at subspot resolution with BayesSapace”. In: Nature biotechnology 39.11 (2021), pp. 1375–1384
- [3] Xiao Zhang et al. “iSC. MEB: an R package for multi-sample spatial clustering analysis of spatial transcriptomics data”. In: Bioinformatics Advances 3.1 (2023), vbad019.
- [4] Carter Allen et al. “MAPLE: a hybrid framework for multi-sample spatial transcriptomics data”. In: bioRxiv (2022), pp. 2022–02.
- [5] Seurat. Analysis, visualization, and integration of spatial datasets with Seurat. https://satijalab.org/seurat/articles/spatial_vignette.html. Accessed: 2023-10-31.
- [6] Simone Avesani et al. “Stardust: improving spatial transcriptomics data analysis through space-aware modularity optimization-based clustering”. In: GigaScience 11 (2022), giac075.
- [7] Yahui Long et al. “Spatially informed clustering, integration, and deconvolution of spatial transcriptomics with GraphST”. In: Nature Communications 14.1 (2023), p. 1155.
- [8] Jian Hu et al. “SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network”. In: Nature methods 18.11 (2021), pp. 1342–1351.
- [9] Kangning Dong and Shihua Zhang. “Deciphering spatial domains from spatially resolved transcriptomics with an adaptive graph attention auto-encoder”. In: Nature communications 13.1 (2022), p. 1739.
- [10] Duy Pham et al. “stLearn: integrating spatial location, tissue morphology and gene expression to find cell types, cell-cell interactions and spatial trajectories within undissociated tissues”. In: BioRxiv (2020), pp. 2020–05.

