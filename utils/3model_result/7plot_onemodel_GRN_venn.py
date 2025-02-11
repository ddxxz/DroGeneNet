import matplotlib.pyplot as plt
import sys
sys.path.append(r'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/src/pyvenn')
import venn
from matplotlib_venn import venn3
import pandas as pd
import string
from icecream import install,ic
install()
import random
config = {
    "font.family":'Times New Roman',
    "font.size": 22,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)
import numpy as  np
data_names=['indica','zeamays','japonica',]
#species=['indica_BGI','japonica','zeamays']
data_types=['TF','PPI','KEGG']#
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']

    for i,speci in enumerate(species):
        # leaf_indica_up_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/leaf_{speci}_CK_CD_DEG_up.csv')
        # leaf_indica_down_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/leaf_{speci}_CK_CD_DEG_down.csv')
        # root_indica_up_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/root_{speci}_CK_CD_DEG_up.csv')
        # root_indica_down_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/root_{speci}_CK_CD_DEG_down.csv')
        GRN_CK_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CK_all.txt',sep='\t')
        GRN_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CD_all.txt',sep='\t')

        #CK_gene_tf=GRN_CK_data['TF']
        #CK_gene_target=GRN_CK_data['Target']
        CK_gene=pd.concat([GRN_CK_data['TF'],GRN_CK_data['Target']],axis=0)

        #CD_gene=GRN_CD_data['FEATURE_ID']
        CD_gene=pd.concat([GRN_CD_data['TF'],GRN_CD_data['Target']],axis=0)
        # 示例数据集
        # 假设三个数据集合为 set1, set2, set3

        plt.figure()

        set1=set(CK_gene)
        set2=set(CD_gene)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn2(labels,figsize=(6,6),fontsize=22, names=[f'{data_names[i]} GRN CK',f'{data_names[i]} GRN CD'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} GRN',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_GRN2.png',bbox_inches='tight')
        plt.clf()


    import numpy as  np
    for i,speci in enumerate(species):
        # leaf_indica_up_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/leaf_{speci}_CK_CD_DEG_up.csv')
        # leaf_indica_down_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/leaf_{speci}_CK_CD_DEG_down.csv')
        # root_indica_up_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/root_{speci}_CK_CD_DEG_up.csv')
        # root_indica_down_data=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_organ_DEG/root_{speci}_CK_CD_DEG_down.csv')
        GRN_CK_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CK_all.txt',sep='\t')
        GRN_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CD_all.txt',sep='\t')
        GRN_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CE_all.txt',sep='\t')

        CK_gene=pd.concat([GRN_CK_data['TF'],GRN_CK_data['Target']],axis=0)
        CD_gene=pd.concat([GRN_CD_data['TF'],GRN_CD_data['Target']],axis=0)
        CE_gene=pd.concat([GRN_CE_data['TF'],GRN_CE_data['Target']],axis=0)
        # 示例数据集
        # 假设三个数据集合为 set1, set2, set3

        plt.figure()

        set1=set(CK_gene)
        set2=set(CD_gene)
        set3=set(CE_gene)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2, set3], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn3(labels,figsize=(6,6),fontsize=22,  names=[f'{data_names[i]} GRN CK',f'{data_names[i]} GRN CD',f'{data_names[i]} GRN CE'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} GRN',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_GRN3.png',bbox_inches='tight')
        plt.clf()













