import numpy as np
import pandas as pd
from icecream import install,ic
install()
import pickle

indica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica/indica_kegg_gene.csv')

indica_tf = indica_regulation.iloc[:,0]#.map(gene_map)
indica_target = indica_regulation.iloc[:,1]#.map(gene_map)
indica_reg=pd.concat([indica_tf,indica_target],axis=1)

japonica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/japonica_kegg_gene.csv')
japonica_tf = japonica_regulation.iloc[:,0]#.apply(convert_msu_to_rap)
japonica_target = japonica_regulation.iloc[:,1]#.apply(convert_msu_to_rap)
japonica_reg=pd.concat([japonica_tf,japonica_target],axis=1)

zeamays_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_kegg_gene.csv')
map_txt2=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo_3/zeamay_geneid_trans.csv',sep='\t')
gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
zeamays_tf2 = zeamays_regulation.iloc[:,0].map(gene_map2)
zeamays_target2 = zeamays_regulation.iloc[:,1].map(gene_map2)
zeamays_reg=pd.concat([zeamays_tf2,zeamays_target2],axis=1)

wheat_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/wheat_kegg_gene.csv')
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
gene_map1 = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
wheat_tf = wheat_regulation.iloc[:,0].map(gene_map1)
wheat_target = wheat_regulation.iloc[:,1].map(gene_map1)
wheat_reg=pd.concat([wheat_tf,wheat_target],axis=1)

sorghum_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/sorghum/sorghum_kegg_gene.csv')

map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.txt')
gene_map1 = dict(zip(map_txt['STRING ID'].astype(str).str.split('.').str[1], map_txt['Gene stable ID']))

map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
gene_map2 = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))

sorghum_tf1 = sorghum_regulation.iloc[:,0].map(gene_map1)
sorghum_target1 = sorghum_regulation.iloc[:,1].map(gene_map1)
sorghum_tf2 = sorghum_tf1.map(gene_map2)
sorghum_target2 = sorghum_target1.map(gene_map2)
sorghum_reg=pd.concat([sorghum_tf2,sorghum_target2],axis=1)
ic(sorghum_reg)
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

indica_reg = indica_reg.drop_duplicates(subset=['TF', 'Target'])
japonica_reg = japonica_reg.drop_duplicates(subset=['TF', 'Target'])
wheat_reg = wheat_reg.drop_duplicates(subset=['TF', 'Target'])
zeamays_reg = zeamays_reg.drop_duplicates(subset=['TF', 'Target'])
sorghum_reg = sorghum_reg.drop_duplicates(subset=['TF', 'Target'])
# ic(zeamays_reg)
# ic(sorghum_reg)
ic(indica_reg)
ic(japonica_reg)
ic(zeamays_reg)
common_reg = pd.merge(indica_reg, japonica_reg, on=['TF', 'Target'], how='inner')
ic(common_reg)
common_reg = pd.merge(common_reg, zeamays_reg, on=['TF', 'Target'], how='inner')
ic(common_reg)
common_reg = pd.merge(common_reg, wheat_reg, on=['TF', 'Target'], how='inner')
ic(common_reg)
common_reg = pd.merge(common_reg, sorghum_reg, on=['TF', 'Target'], how='inner')
ic(common_reg)
# common_reg = pd.merge(zeamays_reg, sorghum_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, indica_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, japonica_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
# common_reg = pd.merge(common_reg, wheat_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)

gene=pd.concat([common_reg['TF'],common_reg['Target']],axis=0)
unique_gene = np.unique(gene)
common_reg.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/homo/TF_Target.csv',index=False)















