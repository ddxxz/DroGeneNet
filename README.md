# DroGeneNet
Training/inference pipeline submission for crop drought-resistance gene regulatory networks **Conserved and divergent gene regulatory networks for crop drought resistance**.

# DroGeneNet Model Framework
The Overall Architecture of the DroGeneNet Model
![图片1](https://github.com/ddxxz/DroGeneNet/blob/main/pic/DroGeneNet.png)
Gene expression embedding layer, large-scale adjacency matrix embedding layer and drought stress encoding layer
![图片2](https://github.com/ddxxz/DroGeneNet/blob/main/pic/DroGeneNet-1.png)
The core component of the graph attention mechanism takes into account the overall connectivity of the network.
![图片3](https://github.com/ddxxz/DroGeneNet/blob/main/pic/DroGeneNet-2.png)
The core component of the graph attention mechanism takes into account the overall connectivity of the network.
Comparison of accuracy for each model in TRN.
![图片4](https://github.com/ddxxz/DroGeneNet/blob/main/pic/model-result.png)
Comparison of accuracy for each model in TRN、PPI、KEGG.
| Layer  | Model      | Japonica AUC | Japonica AUPR | Indica AUC | Indica AUPR | Wheat AUC | Wheat AUPR | Maize AUC | Maize AUPR | Sorghum AUC | Sorghum AUPR |
|--------|------------|--------------|---------------|------------|-------------|-----------|------------|-----------|------------|-------------|--------------|
| TRN    | GENIE3     | 0.50         | 0.50          | 0.50       | 0.50        | 0.50      | 0.49       | 0.50      | 0.49       | 0.52        | 0.52         |
|        | GNE        | 0.92         | 0.92          | 0.91       | 0.93        | 0.85      | 0.82       | 0.88      | 0.80       | 0.83        | 0.82         |
|        | CNNC       | 0.75         | 0.69          | 0.91       | 0.90        | 0.84      | 0.82       | 0.75      | 0.69       | 0.84        | 0.83         |
|        | GENELink   |              |               |            |             |           |            |           |            |             |              |
|        | GNNLink    | 0.75         | 0.69          | 0.75       | 0.69        | 0.82      | 0.86       | 0.74      | 0.68       | 0.76        | 0.66         |
|        | DroGeneNet | 0.94         | 0.93          | 0.94       | 0.93        | 0.89      | 0.88       | 0.86      | 0.84       | 0.85        | 0.83         |
| PPI    | GENIE3     | 0.53         | 0.52          | 0.55       | 0.54        | 0.57      | 0.55       | 0.60      | 0.54       | 0.57        | 0.52         |
|        | GNE        | 0.90         | 0.86          | 0.90       | 0.86        | 0.88      | 0.86       | 0.87      | 0.83       | 0.87        | 0.85         |
|        | CNNC       | 0.76         | 0.71          | 0.75       | 0.69        | 0.78      | 0.70       | 0.79      | 0.77       | 0.76        | 0.72         |
|        | GENELink   | 0.91         | 0.87          | 0.92       | 0.87        | 0.90      | 0.85       | 0.90      | 0.79       | 0.89        | 0.81         |
|        | GNNLink    | 0.75         | 0.69          | 0.75       | 0.69        | 0.76      | 0.68       | 0.75      | 0.69       | 0.75        | 0.69         |
|        | DroGeneNet | 0.93         | 0.88          | 0.93       | 0.88        | 0.97      | 0.97       | 0.90      | 0.84       | 0.94        | 0.89         |
| KEGG   | GENIE3     | 0.52         | 0.50          | 0.54       | 0.52        | 0.54      | 0.52       | 0.56      | 0.54       | 0.54        | 0.50         |
|        | GNE        | 0.85         | 0.88          | 0.88       | 0.90        | 0.88      | 0.87       | 0.77      | 0.79       | 0.82        | 0.77         |
|        | CNNC       | 0.74         | 0.71          | 0.79       | 0.77        | 0.79      | 0.72       | 0.75      | 0.69       | 0.77        | 0.72         |
|        | GENELink   | 0.94         | 0.94          | 0.95       | 0.95        | 0.94      | 0.91       | 0.92      | 0.92       | 0.85        | 0.84         |
|        | GNNLink    | 0.75         | 0.69          | 0.75       | 0.69        | 0.75      | 0.68       | 0.76      | 0.70       | 0.74        | 0.69         |
|        | DroGeneNet | 0.96         | 0.96          | 0.96       | 0.96        | 0.97      | 0.96       | 0.93      | 0.93       | 0.89        | 0.88         |


## Overview
```
├── conf                      <- Data Configuration
│   ├── expression             <- gene expression data
│   ├── link                  <- gene regulatory data
│   ├── config.yaml                  <- Model parameters Configuration
│
├── data                       <-  Project data
│   ├── (ExpressionData_unique_networkfilter.pkl)         <- gene expression data
│   ├── (Train_set.csv)         <- gene regulatory train data
│   ├── (Validation_set.csv)         <- gene regulatory Validation data
│   ├── (Test_set.csv)         <- gene regulatory test data
│   ├── (CK_expression.csv)         <- gene expression data under control condition
│   ├── (CD_expression.csv)         <- gene expression data under drought condition
│
├── data_proceed                 <-  Model data load
│   ├── data_graph.py                 <- Model data load file
├── models                 <-  Model architecture
│   ├── Model_graph.py                 <- Model
├── scripts                        <- Training Run
│   ├── train_all.sh                      <- Training Run file
├── out
│   ├── exp_Graphormer2Link_gate_type_all_type_TF_japonica_japonica                      <- GRN model trained on japonica
│   ├── exp_Graphormer2Link_gate_type_all_type_TF_indica_BGI_indica_BGI                      <- GRN model trained on indica
│   ├── exp_Graphormer2Link_gate_type_all_type_TF_wheat_wheat                      <- GRN model trained on wheat
│   ├── exp_Graphormer2Link_gate_type_all_type_TF_zeamays_zeamays                      <- GRN model trained on zeamays
│   ├── exp_Graphormer2Link_gate_type_all_type_TF_sorghum_sorghum                     <- GRN model trained on sorghum
│
├── utils                 <-  Data Processing and Network Analysis  
│
├── .gitignore                       <- List of files ignored by git
├── environment.yaml                 <- Conda environment file
├── main_graph.py                      <- Main Script
└── README.md
```

File shared via cloud storage: DroGeneNet_data
link: https://pan.baidu.com/s/19f_FuV9LDYi6iz6yUoNGiw?pwd=rwjy code: rwjy 

## How to run

### Clone the repo
```
git clone https://github.com/ddxxz/DroGeneNet.git && cd DroGeneNet
```
### Set up conda environment
This will create conda environment named `DroGeneNet`.
```
mamba env create -n DroGeneNet -f environment.yml
or 
conda env create -n DroGeneNet -f environment.yml
```
### Activate conda environment
```
conda activate DroGeneNet
```
### Training
```
bash train_all.sh
```
### Inference on test data
```
python DroGeneNet/utils/3model_result/2plot_model_weight_onemodel_use.py
python DroGeneNet/utils/3model_result/3GRN_graph_onemodel_use.py
```
### Postanalysis
The specific code can be found in the utils folder

