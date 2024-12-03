## QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport

This repository holds the official source codes of the **QOT** package for the paper [QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport]()

```
@article {QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport
Zexuan Wang, Qipeng Zhan, Shu Yang, Shizhuo Mu, Jiong Chen, Sumita Garai, Patryk Orzechowski, Joost Wagenaar, Li Shen
bioRxiv 2024.02.06.578032; doi: https://doi.org/10.1101/2024.02.06.578032
}
```

## 💡 Abstract
Single-cell technologies have emerged as a transformative technology enabling high-dimensional characterization of cell populations at an unprecedented scale. The data's innate complexity and voluminous nature pose significant computational and analytical challenges, especially in comparative studies delineating cellular architectures across various biological conditions (i.e., generation of sample level distance matrices). Optimal Transport (OT) is a mathematical tool that captures the intrinsic structure of data geometrically and has been applied to many bioinformatics tasks. In this paper, we propose QOT (Quantized Optimal Transport), a new method enables efficient computation of sample level distance matrix from large-scale single-cell omics data through a quantization step. We apply our algorithm to real-world single-cell genomics and pathomics datasets, aiming to extrapolate cell-level insights to inform sample level categorizations. Our empirical study shows that QOT outperforms OT-based algorithms in terms of accuracy and robustness when obtaining a distance matrix at the sample level from high throughput single-cell measures. Moreover, the sample level distance matrix could be used in downstream analysis (i.e. uncover the trajectory of disease progression), highlighting its usage in biomedical informatics and data science.
![alt text](https://github.com/PennShenLab/QOT/blob/main/flow.png)


## 💡 Data
Simulation Datasets could be found at the github folder under Simulation Datasets.
Real-world Datasets could be downloaded at [here](https://zenodo.org/records/8370081) and [here](https://zenodo.org/records/7957118).


## 💡 Requirment
Installation is explicitly handled within each notebook, and all experiments are executed through Google Colab. This setup eliminates the need for users to manage tedious requirements. Additionally, the application of subgroup detection partially relies on the [PILOT](https://github.com/CostaLab/PILOT/tree/main) package, which is also installed directly in the code. 


## 💡 API Reference

The QOT method takes standard [Annotated data](https://anndata.readthedocs.io/en/stable/). It should work well if your datasets follow the same criteria.

### API for Run_QOT

| Parameter             | Description                                                                 |
| :-------------------- | :-------------------------------------------------------------------------- |
| `adata`               | Anndata file for single-cell sequencing                                     |
| `gene_matrix`         | Data matrix                                                                 |
| `type_cell`           | Cell type                                                                   |
| `id_col`              | Sample ID                                                                   |
| `progression`         | Disease status                                                              |
| `dataset_type`        | RNA: use `adata.obsm[gene_matrix]` or Pathomics: use `adata.X`               |
| `num_components_list` | Number of components for each GMM, integer number                           |
| `random_state`        | Random state for reproducibility                                           |
| `min_samples_for_gmm` | Specify minimum samples to consider GMM (can be used as denoise). Default is 1 |
| `qot_method`          | Ground metric for QOT algorithm, choices: "cosine", "euclidean", "exact"   |
| `normalized_set`      | Whether to normalize the QOT distance matrix                               |

### API for trajectory_analysis

| Parameter      | Description                                        |
| :------------- | :------------------------------------------------- |
| `adata`        | Anndata file for single-cell sequencing            |
| `knn`          | Nearest neighbors for building diffusion graph     |
| `dataset_name` | Input dataset to access real labels                |

### API for compute_shapley_values

| Parameter             | Description                                                                 |
| :-------------------- | :-------------------------------------------------------------------------- |
| `adata`               | Anndata file for single-cell sequencing                                     |
| `n_components`        | Number of PCA components for computing Shapley values                       |
| `num_components_list` | Number of components for each GMM, integer number                           |
| `random_state`        | Random state for reproducibility                                           |
| `min_samples_for_gmm` | Specify minimum samples to consider GMM (can be used as denoise). Default is 1 |
| `qot_method`          | Ground metric for QOT algorithm, choices: "cosine", "euclidean", "exact"   |
| `gene_matrix`         | Data matrix                                                                 |
| `type_cell`           | Cell type                                                                   |
| `id_col`              | Sample ID                                                                   |
| `progression`         | Disease status                                                              |
| `dataset_type`        | RNA: use `adata.obsm[gene_matrix]` or Pathomics: use `adata.X`               |



## 💡 Usage
We provide two Jupyter notebooks containing all the experiments and results discussed in our paper. To run each notebook, ensure that the qot_utils.re.py file is placed under the content directory in Google Colab. Additionally, update the file paths for each dataset accordingly:
```terminal
QOT_PDAC_Example.ipynb
```
This notebook first computes the sample-level distance matrix using our QOT algorithm. It then conducts three biological experiments:
- Trajectory Inference
- Driven Gene Identification
- Subgroup Detection, along with differential expression analysis.

```terminal
QOT_Rest.ipynb 
```
This notebook contains the remaining results.

## 💡 Contacts

- [Zexuan Wang](mailto:zxwang@sas.upenn.edu) 
- [Li Shen](mailto:li.shen@pennmedicine.upenn.edu) 

