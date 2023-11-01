# Efficient Mining of Mass Cytometry Data through Optimal Transport

This repository holds the official source codes of the QOT_CYTO (Quantized Optimal Transport for Mass Cytometry) package for the paper [Efficient Mining of Mass Cytometry Data through Optimal Transport]()

### Abstract
Mass Cytometry has emerged as a transformative technology enabling high-dimensional characterization of cell populations at an unprecedented scale. It offers over 40 channels capturing intricate cellular features, aiding in the comprehensive analysis of immune responses and cellular het- erogeneity in health and disease. However, the innate complexity and the volume of data generated present computational and analytical challenges, especially in comparative studies aiming to delineate cellular architectures across different biological conditions. Optimal Transport (OT) is a useful tool to capture the intrinsic structure of data geometrically. Algorithms based on linear programming for calculating OT distances have a cubic time complexity with respect to the input size. This characteristic renders OT inefficient and often infeasible for scenarios involving large samples. In this paper, we present a novel algorithm that enables efficient computation of large-scale mass cytometry data through a quantization step. We apply our algorithm to mass cytometry data derived from COVID-19 patients, aiming to extrapolate cell-level insights to inform population-level categorizations. Our empirical study yields promising results, showing the effectiveness of the proposed algorithm in significantly expedited computation in comparison with multiple competing methods. Additionally, it accurately approximates the original OT-based Wasserstein distance, guaranteeing high fidelity in the representation of intrinsic data characteristics.

### Requirements
QOT_CYTO is implemented in Python. To install the dependent packages, use
```bash
pip install -r requirements.txt
```

### Usage

### Correspondances

[Zexuan Wang](mailto:zxwang@sas.upenn.edu) 
[Li Shen](mailto:li.shen@pennmedicine.upenn.edu) 
