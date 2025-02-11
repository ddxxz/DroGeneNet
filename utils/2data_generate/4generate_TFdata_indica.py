import pandas as pd
import os
import numpy as np
import pickle
from icecream import install,ic
install()

file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_mapos.txt'

# 读取TXT文件，假设以制表符分隔，并没有表头信息
df = pd.read_csv(file_path, usecols=[0, 2])

# 重命名列以方便引用
df.columns = ['Column1', 'Column3']

# 打印原始DataFrame
# print("Original DataFrame:")
# print(df)

# 删除第三列中有缺失值的行
df = df.dropna(subset=['Column3'])

df.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_os.csv',index=False)


import pandas as pd

# 替换为你的映射文件和数据文件路径
map_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/indica_gene_os.csv'
data_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica_BGI/indica_regulation_merged_Osi.txt'

# 读取映射关系文件
map_df = pd.read_csv(map_file_path)
# 创建映射字典，假设映射文件中，第一列为原始基因名，第二列为映射后的基因名
map_dict = pd.Series(map_df.iloc[:, 1].values, index=map_df.iloc[:, 0]).to_dict()

# 读取需要映射的数据文件
data_df = pd.read_csv(data_file_path, sep='\t', header=None)
# 选择第一列和第三列，列索引从0开始
data_df = data_df.iloc[:, [0, 2]]
# 为选中的列设置新的列名
data_df.columns = ['TF', 'Target']

# 假设需要映射的列是第一列和第二列
data_df['Mapped_Column1'] = data_df.iloc[:, 0].map(map_dict)
data_df['Mapped_Column2'] = data_df.iloc[:, 1].map(map_dict)

# 删除任一映射列中含有 NaN 的行
data_df = data_df.dropna(subset=['Mapped_Column1', 'Mapped_Column2'])
data_df = data_df.drop(columns=['TF','Target'])
# 保存或显示结果
print(data_df)
data_df.columns = ['TF','Target']
data_df = data_df.drop_duplicates(subset=['TF', 'Target'])
# 如果需要保存处理后的数据到新的 CSV 文件
data_df.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/TF_Target.csv', index=False)
df_selected = data_df


combined_series = pd.concat([df_selected['TF'], df_selected['Target']])
# 获取唯一值
unique_values = combined_series.unique()
# 转换为DataFrame，为了保存到CSV
unique_df = pd.DataFrame(unique_values, columns=['GeneID'])
# 保存到CSV文件
output_csv_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/regulate_unique_genes.csv'  # 设定输出文件名
unique_df.to_csv(output_csv_path, index=False)


# unique_gene_df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/indica/regulate_unique_genes.csv')
# ic(unique_gene_df.shape)
# #unique_gene_df =unique_gene_df.iloc[:gene_num,]
# # 读取expression data的CSV文件
# with open('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/exp_fpkm_gene_Oryza_indica.pkl', 'rb') as file:
#     # 使用 pickle 加载数据
#     expression_data_df = pickle.load(file)

# # 根据unique gene中的基因符号筛选 expression data
# filtered_expression_data_df = expression_data_df[expression_data_df['FEATURE_ID'].isin(unique_gene_df['GeneID'])]
# cols = filtered_expression_data_df.columns[1:]
# filtered_expression_data_df[cols] = filtered_expression_data_df[cols].applymap(lambda x: x if x >= 0 else 0)
# ic(len(filtered_expression_data_df))
# ic(filtered_expression_data_df)
# ic(filtered_expression_data_df.shape)
# # 保存筛选后的结果为新的CSV文件
# filtered_expression_data_df.to_pickle(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/indica/ExpressionData_unique_networkfilter.pkl')



#-----------------------------------222生成TF索引文件------------------------------------------第一步  与调控关系交叉会减少基因数目
import pandas as pd

df = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/TF_Target.csv")

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/indica_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    Expressiondata = pickle.load(file)
#Expressiondata = Expressiondata.iloc[:11000,:]
Expressiondata.reset_index(inplace=True)
ic(len(df))
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
ic(unique_expressiongenes_sorted)
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
merged_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/ALL_Genefile.csv')



#print(target)
pd.DataFrame({'TF': list(dict.fromkeys(selected_rows['TF'])), 'index': list(dict.fromkeys(tf))}).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/TF.csv')
pd.DataFrame({'Gene': list(dict.fromkeys(selected_rows['Target'])), 'index': list(dict.fromkeys(target))}).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/Target.csv')

pd.DataFrame(unique_genes).to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/unique_genes.csv', index=False)
# 保存转换后的索引对索引的文件
ic(data)
data.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/index_mapped_TF_Target.csv', index=False)

# # # #检查expressiongene得长度


# #---------------------------------------------------555生成label负样本随机采样12万-------------------------------------------------第四步
import pandas as pd
import itertools
import random
from tqdm import tqdm
from icecream import install,ic
install()

# 读取已有的CSV文件
df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/index_mapped_TF_Target.csv')
gene_file = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/ALL_Genefile.csv')
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
negative_samples = random.sample(remaining_pairs, 3*len(existing_pairs))#remaining_pairs ##random.sample(remaining_pairs, 3*len(existing_pairs))
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
all_samples.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/all_labels.csv',index=False)

# # #--------------------------------------444保存uniquegene与网络筛选的表达数据---------------------------------------------------2086 第三步

# import pandas as pd

# 读取unique gene的CSV文件
unique_gene_df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/unique_genes.csv')
#print(unique_gene_df.shape)
#unique_gene_df =unique_gene_df.iloc[:10000,]
# 读取expression data的CSV文件
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/indica_expression_batchcorrection.pkl', 'rb') as file:
    # 使用 pickle 加载数据
    expression_data_df = pickle.load(file)
#expression_data_df = pd.read_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv')
expression_data_df.reset_index(inplace=True)
# 根据unique gene中的基因符号筛选 expression data
filtered_expression_data_df = expression_data_df[expression_data_df['FEATURE_ID'].isin(unique_gene_df['0'])]
# cols = filtered_expression_data_df.columns[1:]
# filtered_expression_data_df[cols] = filtered_expression_data_df[cols].applymap(lambda x: x if x >= 0 else 0)
filtered_expression_data_df = filtered_expression_data_df.drop_duplicates(subset='FEATURE_ID', keep='first')
print(len(filtered_expression_data_df))
# 保存筛选后的结果为新的CSV文件
filtered_expression_data_df.to_pickle('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica/ExpressionData_unique_networkfilter.pkl')
#filtered_expression_data_df.to_csv(f'/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv', index=False)



