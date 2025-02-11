import numpy as np
import pandas as pd
from icecream import install,ic
install()
indica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica_BGI/indica_regulation_merged_Osi.txt',sep='\t')

map_txt=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_mapos.txt')
gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
indica_tf = indica_regulation.iloc[:,0].map(gene_map)
indica_target = indica_regulation.iloc[:,2].map(gene_map)
#ic(indica_tf)
indica_reg=pd.concat([indica_tf,indica_target],axis=1)

japonica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/japonica/japonica_regulation_merged_Osj.txt',sep='\t')
relation={}
for i in open("/home/win/4T/GeneDL/DXZ_DL/expriments/GRN/utils/2data_generate/MSU2RAP/RAP-MSU_2018-03-29.txt"):
    rap=str(i.split()[0])
    msu=str(i.split()[1])
    if msu!="None":
        if "," in msu:
            for a in msu.split(","):
                relation[a[0:-2]] = rap
        else:
            relation[msu[0:-2]] = rap
def convert_rap_to_msu(rap_id):
    return relation.get(rap_id, rap_id)
def convert_msu_to_rap(msu_id):
    return relation.get(msu_id, msu_id)

japonica_tf = japonica_regulation.iloc[:,0].apply(convert_msu_to_rap)
japonica_target = japonica_regulation.iloc[:,2].apply(convert_msu_to_rap)
japonica_reg=pd.concat([japonica_tf,japonica_target],axis=1)


zeamays_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/zeamays/regulation_merged_Zma.txt',sep='\t')
file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/zeamays/gene_id_trans.txt'
map_txt1 = pd.read_csv(
    file_path,
    sep='\t',
    usecols=[0, 1],  # 只读取前两列
    names=['input', '0'],  # 假设文件没有表头，手动指定列名
    header=None  # 告诉pandas没有表头
)
#map_txt1=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt',sep='\t')
gene_map1 = dict(zip(map_txt1['input'], map_txt1['0']))
zeamays_tf1 = zeamays_regulation.iloc[:,0].map(gene_map1)
zeamays_target1 = zeamays_regulation.iloc[:,2].map(gene_map1)

map_txt2=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo_3/zeamay_geneid_trans.csv',sep='\t')
gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
zeamays_tf2 = zeamays_tf1.map(gene_map2)
zeamays_target2 = zeamays_target1.map(gene_map2)
zeamays_reg=pd.concat([zeamays_tf2,zeamays_target2],axis=1)

wheat_regulation=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/regulation_merged_Tae.txt',sep='\t')
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_geneid_trans.csv')
gene_map1 = dict(zip(map_txt['IWGSC1 CSS'], map_txt['IWGSC RefSeq v1.1']))

map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
gene_map2 = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))

wheat_tf1 = wheat_regulation.iloc[:,0].map(gene_map1)
wheat_target1 = wheat_regulation.iloc[:,2].map(gene_map1)
wheat_tf = wheat_tf1.map(gene_map2)
wheat_target = wheat_target1.map(gene_map2)
wheat_reg=pd.concat([wheat_tf,wheat_target],axis=1)
ic(wheat_reg)

sorghum_regulation=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/regulation_merged_Sbi.txt',sep='\t')

map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv')
gene_map1 = dict(zip(map_txt['merged_value'], map_txt['Protein stable ID']))

map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
gene_map2 = dict(zip(map_txt['Protein stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))

sorghum_tf1 = sorghum_regulation.iloc[:,0].map(gene_map1)
sorghum_target1 = sorghum_regulation.iloc[:,2].map(gene_map1)
sorghum_tf = sorghum_tf1.map(gene_map2)
sorghum_target = sorghum_target1.map(gene_map2)
sorghum_reg=pd.concat([sorghum_tf,sorghum_target],axis=1)

#ic(indica_reg,japonica_reg,zeamays_reg)

indica_reg.columns=['TF','Target']
japonica_reg.columns=['TF','Target']
zeamays_reg.columns=['TF','Target']
wheat_reg.columns=['TF','Target']
sorghum_reg.columns=['TF','Target']

zeamays_reg = zeamays_reg.dropna()
sorghum_reg = sorghum_reg.dropna()
wheat_reg = wheat_reg.dropna()
indica_reg = indica_reg.dropna()
japonica_reg = japonica_reg.dropna()

ic(indica_reg)
indica_reg = indica_reg.drop_duplicates(subset=['TF', 'Target'])
ic(indica_reg)
japonica_reg = japonica_reg.drop_duplicates(subset=['TF', 'Target'])
wheat_reg = wheat_reg.drop_duplicates(subset=['TF', 'Target'])
zeamays_reg = zeamays_reg.drop_duplicates(subset=['TF', 'Target'])
sorghum_reg = sorghum_reg.drop_duplicates(subset=['TF', 'Target'])

# ic(zeamays_reg)
# ic(sorghum_reg)
# common_reg = pd.merge(indica_reg, japonica_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, zeamays_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, wheat_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, sorghum_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
common_reg = pd.merge(zeamays_reg, sorghum_reg, on=['TF', 'Target'], how='inner')
#ic(common_reg)
# common_reg = pd.merge(common_reg, indica_reg, on=['TF', 'Target'], how='inner')
# #ic(common_reg)
# common_reg = pd.merge(common_reg, japonica_reg, on=['TF', 'Target'], how='inner')
# #ic(common_reg)
# common_reg = pd.merge(common_reg, wheat_reg, on=['TF', 'Target'], how='inner')
#ic(common_reg)

gene=pd.concat([common_reg['TF'],common_reg['Target']],axis=0)
unique_gene = np.unique(gene)
common_reg.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo_C4/TF_Target.csv',index=False)















