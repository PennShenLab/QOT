# QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport

This repository holds the official source codes of the **QOT** package for the paper [QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport]()

```
@article {QOT: Efficient Computation of Sample Level Distance Matrix from Single-Cell Omics Data through Quantized Optimal Transport
Zexuan Wang, Qipeng Zhan, Shu Yang, Shizhuo Mu, Jiong Chen, Sumita Garai, Patryk Orzechowski, Joost Wagenaar, Li Shen
bioRxiv 2024.02.06.578032; doi: https://doi.org/10.1101/2024.02.06.578032
}
```

### Abstract
Single-cell technologies have emerged as a transformative technology enabling high-dimensional characterization of cell populations at an unprecedented scale. The data's innate complexity and voluminous nature pose significant computational and analytical challenges, especially in comparative studies delineating cellular architectures across various biological conditions (i.e., generation of sample level distance matrices). Optimal Transport (OT) is a mathematical tool that captures the intrinsic structure of data geometrically and has been applied to many bioinformatics tasks. In this paper, we propose QOT (Quantized Optimal Transport), a new method enables efficient computation of sample level distance matrix from large-scale single-cell omics data through a quantization step. We apply our algorithm to real-world single-cell genomics and pathomics datasets, aiming to extrapolate cell-level insights to inform sample level categorizations. Our empirical study shows that QOT outperforms OT-based algorithms in terms of accuracy and robustness when obtaining a distance matrix at the sample level from high throughput single-cell measures. Moreover, the sample level distance matrix could be used in downstream analysis (i.e. uncover the trajectory of disease progression), highlighting its usage in biomedical informatics and data science.

### Data
Simulation Datasets could be found at the github folder under Simulation Datasets.
Real-world Datasets could be downloaded at [here](https://zenodo.org/records/8370081) and [here](https://zenodo.org/records/7957118).



### Usage
The implementation is based on Python. To check each dataset, simply run the notebook under QOT/Tutorial folder.

### Contacts

- [Zexuan Wang](mailto:zxwang@sas.upenn.edu) 
- [Li Shen](mailto:li.shen@pennmedicine.upenn.edu) 

