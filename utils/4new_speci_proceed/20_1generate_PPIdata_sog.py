import pandas as pd
import os
import numpy as np
import pickle
from icecream import install,ic
install()

# file_path1 = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/regulation_merged_Zma.txt'  # 这里替换为你的文本文件路径
# file_path2 = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/Zma_PROTEIN1_list.txt'  # 这里替换为你的文本文件路径

# # 读取文本文件
# df1 = pd.read_csv(file_path1, sep='\t', header=None)  # 如果文件有表头，去掉 `header=None`
# df2 = pd.read_csv(file_path2, sep='\t', header=None)  # 如果文件有表头，去掉 `header=None`

# # # 选择你需要的两列
# # # 假设我们需要第一列和第二列，列索引从0开始
# column1 = df1.iloc[:, 0]
# column2 = df1.iloc[:, 2]

# column3 = df2.iloc[:, 1]
# combined_data = np.concatenate([df1.iloc[:, 0].values, df1.iloc[:, 2].values, df2.iloc[:, 1].values])

# # 提取唯一值
# unique_values = np.unique(combined_data)
# ic(len(unique_values))
# # 保存到CSV文件
# pd.DataFrame(unique_values).to_csv('/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/genes.csv', index=False)  # 第一列的数据保存为 column1.csv

# files3 = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/gene_id_trans.txt'
# #df3 = pd.read_csv(files3, sep='\t', header=None)  # 如果文件有表头，去掉 `header=None`
# # df3 = pd.read_csv(files3, sep='\t', header=None, engine='python')
# # df3[1] = df3[1].str.split().str[0]
# data = []

# # 使用open函数和with语句读取文件
# with open(files3, 'r') as file:
#     for line in file:
#         # 分割每行的数据，只取前两个值（原始基因名和第一个映射基因名）
#         parts = line.strip().split('\t')
#         if len(parts) >= 2:  # 确保有足够的部分读取
#             original_gene = parts[0]
#             mapped_gene = parts[1]
#             # 保存结果
#             data.append([original_gene, mapped_gene])

# # 将列表转换为DataFrame
# df3 = pd.DataFrame(data, columns=['Original Gene', 'Mapped Gene'])


# # 创建映射字典
# gene_id_map = pd.Series(df3['Mapped Gene'].values, index=df3['Original Gene'].values).to_dict()
# #gene_id_map = pd.Series(df3.iloc[:, 1].values, index=df3.iloc[:, 0].values).to_dict()
# #
# # 替换 column1, column2, column3 中的基因名
# column1_mapped = column1.map(gene_id_map).fillna(column1)  # fillna用于保留未找到映射的原始值
# column2_mapped = column2.map(gene_id_map).fillna(column2)
# column3_mapped = column3.map(gene_id_map).fillna(column3)
# ic(column1_mapped)
# ic(column2_mapped)
# ic(column3_mapped)

# regulation = pd.concat([column1_mapped,column2_mapped],axis=1)
# regulation.columns = ['PROTEIN1','Protein2']
# regulation.to_csv('/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/PROTEIN1_Protein2.csv',index=False)

# PROTEIN1_list = pd.concat([column3_mapped, df2.iloc[:, 2]],axis=1)
# PROTEIN1_list.columns = ['PROTEIN1','Family']
# PROTEIN1_list.to_csv('/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/PROTEIN1_family.csv',index=False)

# filtered_df = regulation[regulation['PROTEIN1'].str.startswith('Zm') & regulation['Protein2'].str.startswith('Zm')]
# filtered_df.to_csv('/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/PROTEIN1_Protein2.csv',index=False)

# #-----------------------------------222生成PROTEIN1索引文件------------------------------------------第一步  与调控关系交叉会减少基因数目
import pandas as pd
#蛋白互作得分阈值选择700（默认400, 低150，高700，极高900，越高可信度越强)
#df = pd.read_csv("/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PROTEIN1/zeamays/PROTEIN1_Protein2.csv")
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/sorghum_PPI_gene.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    df = pickle.load(file)
ic(df)
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/sorghum_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    Expressiondata = pickle.load(file)
Expressiondata.reset_index(inplace=True)
ic(Expressiondata)
mapdata = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv')
gene_id_map1 = pd.Series(mapdata['new_column'].values, index=mapdata['Protein stable ID'].values).to_dict()
Expressiondata['mapped_index'] = Expressiondata['FEATURE_ID'].map(gene_id_map1)
# 删除 mapped_index 为 NaN 的行
Expressiondata_cleaned = Expressiondata.dropna(subset=['mapped_index'])
# 将 mapped_index 作为新的索引
Expressiondata_cleaned['FEATURE_ID'] = Expressiondata_cleaned['mapped_index']
Expressiondata_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')

Expressiondata = Expressiondata_cleaned
ic(Expressiondata)
Expressiongene = Expressiondata['FEATURE_ID']
#ic(Expressiongene)
#unique_expressiongenes = list(set(Expressiongene))

# unique_expressiongenes = []
# seen = set()
# for gene in Expressiongene:
#     if gene not in seen:
#         unique_expressiongenes.append(gene)
#         seen.add(gene)
unique_expressiongenes = list(dict.fromkeys(Expressiongene))


ic(len(Expressiongene))
#ic(unique_expressiongenes)
PROTEIN1 = df['protein1']
Protein2 = df['protein2']
ic(len(PROTEIN1))
ic(len(Protein2))
#ic(df['PROTEIN1'].isin(unique_expressiongenes))
selected_rows = df[df['protein1'].isin(unique_expressiongenes) & df['protein2'].isin(unique_expressiongenes)]
#ic(selected_rows)
#unique_protein2 = list(set(selected_rows['Protein2']))


# unique_protein2 = []
# seen = set()
# for gene in selected_rows['Protein2']:
#     if gene not in seen:
#         unique_protein2.append(gene)
#         seen.add(gene)
ic('1')
unique_protein2 = list(dict.fromkeys(selected_rows['protein2']))

#unique_protein1 = list(set(selected_rows['PROTEIN1']))
# unique_protein1 = []
# seen = set()
# for gene in selected_rows['PROTEIN1']:
#     if gene not in seen:
#         unique_protein1.append(gene)
#         seen.add(gene)
ic('2')
unique_protein1 = list(dict.fromkeys(selected_rows['protein1']))


#ic(len(selected_rows))
#ic(unique_protein2)
#ic(unique_protein1)

# 拼接PROTEIN1列和Protein2列为一列
# genes=[]
# for gene in selected_rows['PROTEIN1']:
#     genes.append(gene)

# for gene in selected_rows['Protein2']:
#     genes.append(gene)
# 假设 selected_rows 是一个Pandas DataFrame
ic('3')
genes = [gene for col in ['protein1', 'protein2'] for gene in selected_rows[col]]



#ic(genes)
# # 获取唯一的基因名称以及按顺序生成索引
#unique_genes = list(set(genes))
# unique_genes = []
# seen = set()
# for gene in genes:
#     if gene not in seen:
#         unique_genes.append(gene)
#         seen.add(gene)
ic('4')
unique_genes = list(dict.fromkeys(genes))
ic('5')
#ic(list(unique_genes))#2086
# 创建一个从元素到其在 Expressiongene 中索引的映射字典
# 如果元素不在 Expressiongene 中，设其索引为 -1
index_map = {gene: idx for idx, gene in enumerate(Expressiongene)}
default_index = -1  # 用于不在列表中的元素

# 使用这个映射字典来进行排序，避免在排序时重复查找索引
sorted_unique_genes = sorted(unique_genes, key=lambda x: index_map.get(x, default_index))

#sorted_unique_genes = sorted(unique_genes, key=lambda x: Expressiongene.tolist().index(x) if x in Expressiongene.tolist() else -1)
# sorted_unique_genes = sorted([item for item in unique_genes if item in Expressiongene.tolist()], 
#                      key=lambda x: Expressiongene.index(x))
ic('6')
gene_index = {gene: idx for idx, gene in enumerate(sorted_unique_genes)}
ic(gene_index)
import json
filename = "/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/gene_index.json"
# 使用json.dump()方法将字典保存到JSON文件
with open(filename, 'w', encoding='utf-8') as f:
    json.dump(gene_index, f, ensure_ascii=False, indent=4)
#ic(gene_index)
# # 根据索引映射替换原始数据框的基因名
ic('7')
# 直接在原 DataFrame 上修改，避免赋值给新变量
# 首先，确定需要保留哪些键
data = selected_rows.copy()
data['protein1'] = data['protein1'].map(gene_index).fillna(data['protein1'])
data['protein2'] = data['protein2'].map(gene_index).fillna(data['protein2'])
#data = selected_rows.replace(gene_index)
ic('8')
#print(gene_index)
protein1 = selected_rows['protein1'].map(gene_index).fillna(selected_rows['protein1'])

#protein1 = selected_rows['PROTEIN1'].replace(gene_index)
#print(protein1)
ic('9')
protein2 = selected_rows['protein2'].map(gene_index).fillna(selected_rows['protein2'])

#protein2 = selected_rows['Protein2'].replace(gene_index)
ic('10')
unique_expressiongenes_sorted = Expressiondata['FEATURE_ID'].unique().tolist()
ic(len(unique_expressiongenes_sorted))
ic('11')
expressiongene_index = {gene: idx for idx, gene in enumerate(unique_expressiongenes_sorted)}
ic('12')
merged_df = pd.concat([pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['protein1'])), 'index': list(dict.fromkeys(protein1))}), pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['protein2'])), 'index': list(dict.fromkeys(protein2))})], axis=0)
merged_df = pd.DataFrame(list(gene_index.items()), columns=['Gene', 'index'])
# 确保merged_df的Gene列按照unique_expressiongenes的顺序排序
# merged_df['sort_index'] = merged_df['Gene'].map(expressiongene_index)
# merged_df_sorted = merged_df.sort_values(by='sort_index').drop('sort_index', axis=1)
#merged_df.columns = ['Gene', 'index']
ic(len(merged_df))
merged_df.to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/ALL_Genefile.csv')



#print(protein2)
pd.DataFrame({'protein1': list(dict.fromkeys(selected_rows['protein1'])), 'index': list(dict.fromkeys(protein1))}).to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/protein1.csv')
pd.DataFrame({'protein2': list(dict.fromkeys(selected_rows['protein2'])), 'index': list(dict.fromkeys(protein2))}).to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/protein2.csv')

pd.DataFrame(unique_genes).to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/unique_genes.csv', index=False)
# 保存转换后的索引对索引的文件

#ic(data['PROTEIN1'].max(),data['Protein2'].max())
data.to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/index_mapped_protein1_protein2.csv', index=False)

# # # #检查expressiongene得长度


# #---------------------------------------------------555生成label负样本随机采样12万-------------------------------------------------第四步
import pandas as pd
import itertools
import random
from tqdm import tqdm
from icecream import install,ic
install()
# 读取已有的CSV文件
df = pd.read_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/index_mapped_protein1_protein2.csv')
gene_file = pd.read_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/ALL_Genefile.csv')
df = df.sort_values('protein1', ascending=True)
protein1_unique_values = df['protein1'].unique()
all_combinations =[]

for protein1 in protein1_unique_values:
    combination = [(protein1, i) for i in range(len(gene_file)) if i!=protein1]
    all_combinations.extend(combination)
# 生成无向组合

all_combinations =set(all_combinations)
#ic(df)
df['Label'] = 1
#print(df['PROTEIN1'].max(),df['Protein2'].max())
existing_pairs = df[df['Label'] == 1][['protein1', 'protein2']]
# 生成剩下的排列组合
existing_index_pairs = set(zip(existing_pairs['protein1'], existing_pairs['protein2']))
#all_combinations = set(itertools.combinations(range(10), 2))
#all_combinations = [()]
#ic(all_combinations)
remaining_pairs = all_combinations - existing_index_pairs
ic(len(remaining_pairs))
ic(len(existing_pairs))
#ic(len(remaining_pairs))
#ic(all_combinations)
# 从剩下的索引对中随机抽取负样本---------------------可能有重复
negative_samples = random.sample(remaining_pairs, 3*len(existing_pairs)) #remaining_pairs#
negative_samples = pd.DataFrame(list(set(negative_samples)), columns=['protein1', 'protein2'])
negative_samples = negative_samples.sort_values('protein1', ascending=True)

# 构建正样本和负样本的DataFrame
positive_samples = existing_pairs.copy()
positive_samples['Label'] = 1

negative_samples_df = pd.DataFrame(negative_samples, columns=['protein1', 'protein2'])
negative_samples_df['Label'] = 0

# 合并正样本和负样本，并抽取12万个样本
all_samples = pd.concat([positive_samples, negative_samples_df], ignore_index=True)

# 保存为新的CSV文件
all_samples.to_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/all_labels.csv',index=False)

# # #--------------------------------------444保存uniquegene与网络筛选的表达数据---------------------------------------------------2086 第三步

# import pandas as pd

# 读取unique gene的CSV文件
unique_gene_df = pd.read_csv(f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/unique_genes.csv')
#print(unique_gene_df.shape)
#unique_gene_df =unique_gene_df.iloc[:10000,]
# 读取expression data的CSV文件
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/sorghum_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    expression_data_df = pickle.load(file)
#expression_data_df = pd.read_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv')
expression_data_df.reset_index(inplace=True)
ic(expression_data_df)
mapdata = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv')
gene_id_map1 = pd.Series(mapdata['new_column'].values, index=mapdata['Protein stable ID'].values).to_dict()
expression_data_df['mapped_index'] = expression_data_df['FEATURE_ID'].map(gene_id_map1)
# 删除 mapped_index 为 NaN 的行
Expressiondata_cleaned = expression_data_df.dropna(subset=['mapped_index'])
# 将 mapped_index 作为新的索引
Expressiondata_cleaned['FEATURE_ID'] = Expressiondata_cleaned['mapped_index']
Expressiondata_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')

expression_data_df = Expressiondata_cleaned
# 根据unique gene中的基因符号筛选 expression data
filtered_expression_data_df = expression_data_df[expression_data_df['FEATURE_ID'].isin(unique_gene_df['0'])]
# cols = filtered_expression_data_df.columns[1:]
# filtered_expression_data_df[cols] = filtered_expression_data_df[cols].applymap(lambda x: x if x >= 0 else 0)

print(len(filtered_expression_data_df))
# 保存筛选后的结果为新的CSV文件
filtered_expression_data_df.to_pickle('/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/ExpressionData_unique_networkfilter.pkl')
#filtered_expression_data_df.to_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv', index=False)



