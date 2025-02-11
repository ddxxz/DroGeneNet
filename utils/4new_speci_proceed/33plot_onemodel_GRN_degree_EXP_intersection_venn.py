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

data_names=['wheat','sorghum','homo_C3','homo_C4','homo']
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type == 'TF':
        species=['wheat','sorghum','homo_C3','homo_C4'] 
    else:
        species = ['wheat','sorghum','homo_C3','homo_C4','homo']#
    for i,speci in enumerate(species):

        degree_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CD/all_{speci}_CK_CD_DEG.csv')
        exp_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG.csv')

        degree_CK_CD_deg=degree_CK_CD_data['Unnamed: 0']
        exp_CK_CD_deg=exp_CK_CD_data['FEATURE_ID']

        plt.figure()

        set1=set(degree_CK_CD_deg)
        set2=set(exp_CK_CD_deg)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn2(labels,figsize=(6,6),fontsize=22, names=[f'{data_names[i]} degree CK CD',f'{data_names[i]} exp CK CD'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} CK CD degree&exp',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_CK_CD_degree_exp_DEG.png',bbox_inches='tight')
        plt.clf()

        ck_cd_intersection= set1 & set2
        pd.DataFrame(ck_cd_intersection).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/{data_type}_{speci}_ck_cd_degree_exp_intersection.csv',index=False)

    import numpy as  np
    for i,speci in enumerate(species):

        degree_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CE/all_{speci}_CK_CE_DEG.csv')
        exp_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG.csv')

        degree_CK_CE_deg=degree_CK_CE_data['Unnamed: 0']
        exp_CK_CE_deg=exp_CK_CE_data['FEATURE_ID']

        plt.figure()

        set1=set(degree_CK_CE_deg)
        set2=set(exp_CK_CE_deg)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn2(labels,figsize=(6,6),fontsize=22,names=[f'{data_names[i]} degree CK CE',f'{data_names[i]} exp CK CE'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} CK CE degree&exp',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_CK_CE_degree_exp_DEG.png',bbox_inches='tight')
        plt.clf()

        ck_ce_intersection= set1 & set2
        pd.DataFrame(ck_ce_intersection).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/{data_type}_{speci}_ck_ce_degree_exp_intersection.csv',index=False)

    for i,speci in enumerate(species):

        degree_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CD/all_{speci}_CK_CD_DEG.csv')
        exp_CK_CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG.csv')

        degree_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CE/all_{speci}_CK_CE_DEG.csv')
        exp_CK_CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG.csv')

        degree_CK_CD_deg=degree_CK_CD_data['Unnamed: 0']
        exp_CK_CD_deg=exp_CK_CD_data['FEATURE_ID']

        degree_CK_CE_deg=degree_CK_CE_data['Unnamed: 0']
        exp_CK_CE_deg=exp_CK_CE_data['FEATURE_ID']

        plt.figure()

        set1=set(degree_CK_CD_deg)
        set2=set(exp_CK_CD_deg)
        set3=set(degree_CK_CE_deg)
        set4=set(exp_CK_CE_deg)

        colors = ['#8B0000', '#006400', '#00008B']

        labels = venn.get_labels([set1, set2,set3,set4], fill=['number',
                                                            #'logic',,set3,set4
                                                        #'percent'
                                                                ]
                                )
        fig, ax = venn.venn4(labels,figsize=(6,6),fontsize=22, names=[f'{data_names[i]} degree CK CD',f'{data_names[i]} exp CK CD',f'{data_names[i]} degree CK CE',f'{data_names[i]} exp CK CE'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

        plt.title(f'{data_names[i]} degree&exp',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_venn/{data_type}_{speci}_CK_CD_CE_degree_exp_DEG.png',bbox_inches='tight')
        plt.clf()

        ck_cd_intersection= set1 & set2& set3& set4
        pd.DataFrame(ck_cd_intersection).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/{data_type}_{speci}_ck_cd_ce_degree_exp_intersection.csv',index=False)













