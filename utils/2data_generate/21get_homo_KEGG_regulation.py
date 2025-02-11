import numpy as np
import pandas as pd
from icecream import install,ic
install()
import pickle

indica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica/indica_kegg_gene.csv')
# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/indica_BGI_PPI_gene.pkl','rb') as fp:
#     indica_regulation = pickle.load(fp)
# ic(indica_regulation)
# map_txt=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_mapos.txt')
# gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
indica_tf = indica_regulation.iloc[:,0]#.map(gene_map)
indica_target = indica_regulation.iloc[:,1]#.map(gene_map)
#ic(indica_tf)
indica_reg=pd.concat([indica_tf,indica_target],axis=1)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica_PPI_gene.pkl','rb') as fp:
#     japonica_regulation = pickle.load(fp)

japonica_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/japonica_kegg_gene.csv')
# relation={}
# for i in open("/home/win/4T/GeneDL/DXZ_DL/expriments/GRN/utils/data_generate/MSU2RAP/RAP-MSU_2018-03-29.txt"):
#     rap=str(i.split()[0])
#     msu=str(i.split()[1])
#     if msu!="None":
#         if "," in msu:
#             for a in msu.split(","):
#                 relation[a[0:-2]] = rap
#         else:
#             relation[msu[0:-2]] = rap
# def convert_rap_to_msu(rap_id):
#     return relation.get(rap_id, rap_id)
# def convert_msu_to_rap(msu_id):
#     return relation.get(msu_id, msu_id)

japonica_tf = japonica_regulation.iloc[:,0]#.apply(convert_msu_to_rap)
japonica_target = japonica_regulation.iloc[:,1]#.apply(convert_msu_to_rap)
japonica_reg=pd.concat([japonica_tf,japonica_target],axis=1)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/zeamays_PPI_gene.pkl','rb') as fp:
#     zeamays_regulation = pickle.load(fp)
#ic(zeamays_regulation)
zeamays_regulation=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_kegg_gene.csv')
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/zeamays/gene_id_trans.txt'
# map_txt1 = pd.read_csv(
#     file_path,
#     sep='\t',
#     usecols=[0, 1],  # 只读取前两列
#     names=['input', '0'],  # 假设文件没有表头，手动指定列名
#     header=None  # 告诉pandas没有表头
# )
# #map_txt1=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt',sep='\t')
# gene_map1 = dict(zip(map_txt1['input'], map_txt1['0']))
# zeamays_tf1 = zeamays_regulation.iloc[:,0].map(gene_map1)
# zeamays_target1 = zeamays_regulation.iloc[:,1].map(gene_map1)
#ic(gene_map1)
map_txt2=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo/zeamay_geneid_trans.csv',sep='\t')
gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
zeamays_tf2 = zeamays_regulation.iloc[:,0].map(gene_map2)
zeamays_target2 = zeamays_regulation.iloc[:,1].map(gene_map2)
zeamays_reg=pd.concat([zeamays_tf2,zeamays_target2],axis=1)

#ic(indica_reg,japonica_reg,zeamays_reg)
ic(indica_reg)
ic(japonica_reg)
ic(zeamays_reg)
indica_reg.columns=['TF','Target']
japonica_reg.columns=['TF','Target']
zeamays_reg.columns=['TF','Target']


common_reg = pd.merge(indica_reg, japonica_reg, on=['TF', 'Target'], how='inner')

# 将合并结果与 zeamays_reg 进一步合并
common_reg = pd.merge(common_reg, zeamays_reg, on=['TF', 'Target'], how='inner')
gene=pd.concat([common_reg['TF'],common_reg['Target']],axis=0)
unique_gene = np.unique(gene)
ic(common_reg)
ic(len(unique_gene))
common_reg.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/homo/TF_Target.csv',index=False)















