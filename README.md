# WGCNA_CNN
WGCNA_CNN所有的代码
0-Breast_data中包含了针对于TCGA数据的处理，将其处理为基因表达矩阵的形式。对数变换，离群样本剔除也在其中。TCGA数据过大，可以去官网(https://portal.gdc.cancer.gov/)进行下载,在其中勾选breast cancer，然后利用当前代码处理。
1-Screen gene module and key genes中GSE25066文件夹包含了GSE25066的数据获取以及针对于GSE25066数据进行WGCNA+GO的代码；WGCNA_analysis2中包含了针对于TCGA乳腺癌数据进行WGCNA+GO的代码
2-Model中包含了GSE25066_CNN_model_runing（WGCNA+CNN在GSE25066数据上的运行）、Mutual Information（互信息在tcga数据上的运行）、RandomForestClassifier（随机森林在tcga数据上的运行）、yellowModule_WGCNA_CNN（WGCNA+CNN在TCGA数据上的运行）
