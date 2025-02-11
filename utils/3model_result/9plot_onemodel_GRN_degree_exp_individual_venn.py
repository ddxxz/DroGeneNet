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

data_names=['indica','zeamays','japonica']
data_types=['TF','PPI','KEGG']#
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']
    import numpy as  np
    for i,speci in enumerate(species):

        degree_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CD/all_{speci}_CK_CD_DEG.csv')
        degree_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CE/all_{speci}_CK_CE_DEG.csv')
        degree_CD_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CD_CE/all_{speci}_CD_CE_DEG.csv')

        CK_CD_deg=degree_CK_CD_data['Unnamed: 0']
        CK_CE_deg=degree_CK_CE_data['Unnamed: 0']
        CD_CE_deg=degree_CD_CE_data['Unnamed: 0']

        plt.figure()

        set1=set(CK_CD_deg)
        set2=set(CK_CE_deg)
        set3=set(CD_CE_deg)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2,set3], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn3(labels,figsize=(6,6),fontsize=22, names=[f'{data_names[i]} CK CD',f'{data_names[i]} CK CE',f'{data_names[i]} CD CE'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} degree',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_degree_DEG.png',bbox_inches='tight')
        plt.clf()


    import numpy as  np
    for i,speci in enumerate(species):

        degree_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG.csv')
        degree_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG.csv')
        degree_CD_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CD_CE/{speci}_CD_CE_DEG.csv')

        CK_CD_deg=degree_CK_CD_data['FEATURE_ID']
        CK_CE_deg=degree_CK_CE_data['FEATURE_ID']
        CD_CE_deg=degree_CD_CE_data['FEATURE_ID']

        plt.figure()

        set1=set(CK_CD_deg)
        set2=set(CK_CE_deg)
        set3=set(CD_CE_deg)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2,set3], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn3(labels,figsize=(6,6),fontsize=22, names=[f'{data_names[i]} CK CD',f'{data_names[i]} CK CE',f'{data_names[i]} CD CE'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} exp',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_exp_DEG.png',bbox_inches='tight')
        plt.clf()











