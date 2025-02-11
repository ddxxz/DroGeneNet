import numpy as np
import pandas as pd
from icecream import install,ic
install()
import pickle

#indica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica_BGI/indica_regulation_merged_Osi.txt',sep='\t')
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/indica_BGI_PPI_gene.pkl','rb') as fp:
    indica_regulation = pickle.load(fp)
ic(indica_regulation)
map_txt=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_mapos.txt')
gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
indica_tf = indica_regulation.iloc[:,0].map(gene_map)
indica_target = indica_regulation.iloc[:,1].map(gene_map)
#ic(indica_tf)
indica_reg=pd.concat([indica_tf,indica_target],axis=1)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica_PPI_gene.pkl','rb') as fp:
    japonica_regulation = pickle.load(fp)

japonica_tf = japonica_regulation.iloc[:,0]#.apply(convert_msu_to_rap)
japonica_target = japonica_regulation.iloc[:,1]#.apply(convert_msu_to_rap)
japonica_reg=pd.concat([japonica_tf,japonica_target],axis=1)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/zeamays_PPI_gene.pkl','rb') as fp:
    zeamays_regulation = pickle.load(fp)
ic(zeamays_regulation)
map_txt2=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo_3/zeamay_geneid_trans.csv',sep='\t')
gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
zeamays_tf2 = zeamays_regulation.iloc[:,0].map(gene_map2)
zeamays_target2 = zeamays_regulation.iloc[:,1].map(gene_map2)
zeamays_reg=pd.concat([zeamays_tf2,zeamays_target2],axis=1)

wheat_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/wheat/PPI_reg.csv')
map_txt=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/wheat/idmapping_wheat.tsv',sep='\t')
gene_map1 = dict(zip(map_txt['From'], map_txt['To']))
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
gene_map2 = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))

wheat_tf1 = wheat_regulation.iloc[:,0].map(gene_map1)
wheat_target1 = wheat_regulation.iloc[:,1].map(gene_map1)
wheat_tf = wheat_tf1.map(gene_map2)
wheat_target = wheat_target1.map(gene_map2)
wheat_reg=pd.concat([wheat_tf,wheat_target],axis=1)
ic(wheat_reg)

sorghum_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/PPI_reg.csv')
map_txt=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/idmapping_sorghum.tsv',sep='\t')
gene_map1 = dict(zip(map_txt['From'], map_txt['To']))
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
gene_map2 = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))

sorghum_tf1 = sorghum_regulation.iloc[:,0].map(gene_map1)
sorghum_target1 = sorghum_regulation.iloc[:,1].map(gene_map1)
sorghum_tf = sorghum_tf1.map(gene_map2)
sorghum_target = sorghum_target1.map(gene_map2)
sorghum_reg=pd.concat([sorghum_tf,sorghum_target],axis=1)
ic(sorghum_reg)
#ic(indica_reg,japonica_reg,zeamays_reg)
# ic(indica_reg)
# ic(japonica_reg)
# ic(zeamays_reg)
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
#ic(indica_reg)
indica_reg = indica_reg.drop_duplicates(subset=['TF', 'Target'])
#ic(indica_reg)
japonica_reg = japonica_reg.drop_duplicates(subset=['TF', 'Target'])
wheat_reg = wheat_reg.drop_duplicates(subset=['TF', 'Target'])
zeamays_reg = zeamays_reg.drop_duplicates(subset=['TF', 'Target'])
sorghum_reg = sorghum_reg.drop_duplicates(subset=['TF', 'Target'])
# ic(zeamays_reg)
# ic(sorghum_reg)
common_reg = pd.merge(indica_reg, japonica_reg, on=['TF', 'Target'], how='inner')
#ic(common_reg)
# common_reg = pd.merge(common_reg, zeamays_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
common_reg = pd.merge(common_reg, wheat_reg, on=['TF', 'Target'], how='inner')
#ic(common_reg)
# common_reg = pd.merge(common_reg, sorghum_reg, on=['TF', 'Target'], how='inner')
# ic(common_reg)
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
common_reg.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/homo_C3/TF_Target.csv',index=False)














