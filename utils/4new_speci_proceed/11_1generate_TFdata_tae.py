import pandas as pd
import os
import numpy as np
import pickle
from icecream import install,ic
install()

file_path1 = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/regulation_merged_Tae.txt'  # 这里替换为你的文本文件路径
file_path2 = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/Tae_TF_list.txt'  # 这里替换为你的文本文件路径

# 读取文本文件
df1 = pd.read_csv(file_path1, sep='\t', header=None)  # 如果文件有表头，去掉 `header=None`
df2 = pd.read_csv(file_path2, sep='\t', header=None)  # 如果文件有表头，去掉 `header=None`

# # 选择你需要的两列
# # 假设我们需要第一列和第二列，列索引从0开始
column1 = df1.iloc[:, 0]
column2 = df1.iloc[:, 2]
column3 = df2.iloc[:, 1]
combined_data = np.concatenate([df1.iloc[:, 0].values, df1.iloc[:, 2].values, df2.iloc[:, 1].values])

# 提取唯一值
unique_values = np.unique(combined_data)
ic(len(unique_values))
# 保存到CSV文件
pd.DataFrame(unique_values).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/genes.csv', index=False)  # 第一列的数据保存为 column1.csv

files3 = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_geneid_trans.csv'
df3 = pd.read_csv(files3)
data = []

# 创建映射字典
gene_id_map = pd.Series(df3['IWGSC RefSeq v1.1'].values, index=df3['IWGSC1 CSS'].values).to_dict()
#gene_id_map = pd.Series(df3.iloc[:, 1].values, index=df3.iloc[:, 0].values).to_dict()
#
# 替换 column1, column2, column3 中的基因名
column1_mapped = column1.map(gene_id_map).fillna(column1)  # fillna用于保留未找到映射的原始值
column2_mapped = column2.map(gene_id_map).fillna(column2)
column3_mapped = column3.map(gene_id_map).fillna(column3)
ic(column1_mapped)
ic(column2_mapped)
ic(column3_mapped)

regulation = pd.concat([column1_mapped,column2_mapped],axis=1)
regulation.columns = ['TF','Target']
regulation = regulation.drop_duplicates(subset=['TF', 'Target'])
regulation.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/TF_Target.csv',index=False)

TF_list = pd.concat([column3_mapped, df2.iloc[:, 2]],axis=1)
TF_list.columns = ['TF','Family']
TF_list.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/TF_family.csv',index=False)

# filtered_df = regulation[regulation['TF'].str.startswith('Zm') & regulation['Target'].str.startswith('Zm')]
# filtered_df.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/TF_Target.csv',index=False)

# #-----------------------------------222生成TF索引文件------------------------------------------第一步  与调控关系交叉会减少基因数目
import pandas as pd

df = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/TF_Target.csv")

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/wheat_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    Expressiondata = pickle.load(file)
Expressiondata.reset_index(inplace=True)
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
TF = df['TF']
Target = df['Target']
ic(len(TF))
ic(len(Target))
#ic(df['TF'].isin(unique_expressiongenes))
selected_rows = df[df['TF'].isin(unique_expressiongenes) & df['Target'].isin(unique_expressiongenes)]
#ic(selected_rows)
#unique_target = list(set(selected_rows['Target']))


# unique_target = []
# seen = set()
# for gene in selected_rows['Target']:
#     if gene not in seen:
#         unique_target.append(gene)
#         seen.add(gene)
ic('1')
unique_target = list(dict.fromkeys(selected_rows['Target']))

#unique_tf = list(set(selected_rows['TF']))
# unique_tf = []
# seen = set()
# for gene in selected_rows['TF']:
#     if gene not in seen:
#         unique_tf.append(gene)
#         seen.add(gene)
ic('2')
unique_tf = list(dict.fromkeys(selected_rows['TF']))


#ic(len(selected_rows))
#ic(unique_target)
#ic(unique_tf)

# 拼接TF列和Target列为一列
# genes=[]
# for gene in selected_rows['TF']:
#     genes.append(gene)

# for gene in selected_rows['Target']:
#     genes.append(gene)
# 假设 selected_rows 是一个Pandas DataFrame
ic('3')
genes = [gene for col in ['TF', 'Target'] for gene in selected_rows[col]]



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
filename = "/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/gene_index.json"
# 使用json.dump()方法将字典保存到JSON文件
with open(filename, 'w', encoding='utf-8') as f:
    json.dump(gene_index, f, ensure_ascii=False, indent=4)
#ic(gene_index)
# # 根据索引映射替换原始数据框的基因名
ic('7')
# 直接在原 DataFrame 上修改，避免赋值给新变量
# 首先，确定需要保留哪些键
data = selected_rows.copy()
data['TF'] = data['TF'].map(gene_index).fillna(data['TF'])
data['Target'] = data['Target'].map(gene_index).fillna(data['Target'])
#data = selected_rows.replace(gene_index)
ic('8')
#print(gene_index)
tf = selected_rows['TF'].map(gene_index).fillna(selected_rows['TF'])

#tf = selected_rows['TF'].replace(gene_index)
#print(tf)
ic('9')
target = selected_rows['Target'].map(gene_index).fillna(selected_rows['Target'])

#target = selected_rows['Target'].replace(gene_index)
ic('10')
unique_expressiongenes_sorted = Expressiondata['FEATURE_ID'].unique().tolist()
ic(len(unique_expressiongenes_sorted))
ic('11')
expressiongene_index = {gene: idx for idx, gene in enumerate(unique_expressiongenes_sorted)}
ic('12')
merged_df = pd.concat([pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['TF'])), 'index': list(dict.fromkeys(tf))}), pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['Target'])), 'index': list(dict.fromkeys(target))})], axis=0)
merged_df = pd.DataFrame(list(gene_index.items()), columns=['Gene', 'index'])
# 确保merged_df的Gene列按照unique_expressiongenes的顺序排序
# merged_df['sort_index'] = merged_df['Gene'].map(expressiongene_index)
# merged_df_sorted = merged_df.sort_values(by='sort_index').drop('sort_index', axis=1)
#merged_df.columns = ['Gene', 'index']
ic(len(merged_df))
merged_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/ALL_Genefile.csv')



#print(target)
pd.DataFrame({'TF': list(dict.fromkeys(selected_rows['TF'])), 'index': list(dict.fromkeys(tf))}).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/TF.csv')
pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['Target'])), 'index': list(dict.fromkeys(target))}).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/Target.csv')

pd.DataFrame(unique_genes).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/unique_genes.csv', index=False)
# 保存转换后的索引对索引的文件

ic(data['TF'].max(),data['Target'].max())
data.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/index_mapped_TF_Target.csv', index=False)

# # #检查expressiongene得长度


# #---------------------------------------------------555生成label负样本随机采样12万-------------------------------------------------第四步
import pandas as pd
import itertools
import random
from tqdm import tqdm
from icecream import install,ic
install()
# 读取已有的CSV文件
df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/index_mapped_TF_Target.csv')
gene_file = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/ALL_Genefile.csv')
df = df.sort_values('TF', ascending=True)
tf_unique_values = df['TF'].unique()
all_combinations =[]
for tf in tf_unique_values:
    combination = [(tf, i) for i in range(len(gene_file)) if i!=tf]
    all_combinations.extend(combination)
all_combinations =set(all_combinations)
#ic(df)
df['Label'] = 1
#print(df['TF'].max(),df['Target'].max())
existing_pairs = df[df['Label'] == 1][['TF', 'Target']]
# 生成剩下的排列组合
existing_index_pairs = set(zip(existing_pairs['TF'], existing_pairs['Target']))
#all_combinations = set(itertools.combinations(range(10), 2))
#all_combinations = [()]
#ic(all_combinations)
remaining_pairs = all_combinations - existing_index_pairs
ic(len(remaining_pairs))
ic(len(existing_pairs))
#ic(len(remaining_pairs))
#ic(all_combinations)
# 从剩下的索引对中随机抽取负样本---------------------可能有重复
negative_samples = random.sample(remaining_pairs, 3*len(existing_pairs)) #remaining_pairs##
negative_samples = pd.DataFrame(list(set(negative_samples)), columns=['TF', 'Target'])
negative_samples = negative_samples.sort_values('TF', ascending=True)

# 构建正样本和负样本的DataFrame
positive_samples = existing_pairs.copy()
positive_samples['Label'] = 1

negative_samples_df = pd.DataFrame(negative_samples, columns=['TF', 'Target'])
negative_samples_df['Label'] = 0

# 合并正样本和负样本，并抽取12万个样本
all_samples = pd.concat([positive_samples, negative_samples_df], ignore_index=True)

# 保存为新的CSV文件
all_samples.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/all_labels.csv',index=False)

# # #--------------------------------------444保存uniquegene与网络筛选的表达数据---------------------------------------------------2086 第三步
import pickle
# import pandas as pd

# 读取unique gene的CSV文件
unique_gene_df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/unique_genes.csv')
#print(unique_gene_df.shape)
#unique_gene_df =unique_gene_df.iloc[:10000,]
# 读取expression data的CSV文件
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/wheat_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    expression_data_df = pickle.load(file)
#expression_data_df = pd.read_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv')
expression_data_df.reset_index(inplace=True)
# 根据unique gene中的基因符号筛选 expression data
filtered_expression_data_df = expression_data_df[expression_data_df['FEATURE_ID'].isin(unique_gene_df['0'])]
# cols = filtered_expression_data_df.columns[1:]
# filtered_expression_data_df[cols] = filtered_expression_data_df[cols].applymap(lambda x: x if x >= 0 else 0)

print(len(filtered_expression_data_df))
# 保存筛选后的结果为新的CSV文件
filtered_expression_data_df.to_pickle('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/wheat/ExpressionData_unique_networkfilter.pkl')
#filtered_expression_data_df.to_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv', index=False)



