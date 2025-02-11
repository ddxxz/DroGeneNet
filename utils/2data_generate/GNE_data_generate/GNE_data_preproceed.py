import pandas as pd
import numpy as np
import sys
import pickle
# 读取原始 CSV 文件
data_path=sys.argv[1]
species=sys.argv[2]
df = pd.read_csv(f'{data_path}/{species}/ALL_Genefile.csv')
# 创建新的 DataFrame
new_df = pd.DataFrame()
# 将 index 列保存到新的 I#ID 列，并加上字母 G
new_df['I#ID'] = 'G' + (df['index']).astype(str)
# 将 Gene 列保存到 Name 列
new_df['Name'] = df['Gene']
# 打印新的 DataFrame
print(new_df)
# 将 DataFrame 保存到 CSV 文件
new_df.to_csv(f'{data_path}/{species}/gene_ids.tsv', index=False)


# 读取 gene_ids.tsv 文件
gene_ids_df = pd.read_csv(f'{data_path}/{species}/gene_ids.tsv', sep='\s*,\s*', engine='python')
print(gene_ids_df)
# 读取原始 CSV 文件
with open(f'{data_path}/{species}/ExpressionData_unique_networkfilter.pkl', "rb") as fp:
    df = pickle.load(fp)
# 将 Gene Symbol 替换为对应的 I#ID
df = df.dropna(axis=1, how='any')
df['FEATURE_ID'] = df['FEATURE_ID'].map(gene_ids_df.set_index('Name')['I#ID'])
#print(df)
# 转置 dataframe
transposed_df = df.T
# 将 0 替换为 NaN
# # 将非数值列替换为 NaN
# sub_df = transposed_df.iloc[1:]

# cols_with_zero = sub_df.columns[sub_df.isin([0]).any()]
# sub_df[cols_with_zero] = sub_df[cols_with_zero].replace(0, np.nan)

# # 使用上下数的平均值填充 NaN
# sub_df.fillna((sub_df.shift() + sub_df.shift(-1)) / 2, inplace=True)
# transposed_df = pd.concat([transposed_df.iloc[:1], sub_df])
# #transposed_df = transposed_df.iloc[:,1:]
# # transposed_df = pd.DataFrame(columns=transposed_df[0])
print(transposed_df)

# # 保存转置后的 dataframe 到 expression_data.tsv
# transposed_df.to_csv('/mnt/h/study/deep_learning/gene/project/GNE/data/oryza_sativa/expression_data.tsv', sep='\t', header=False, index=False)

transposed_df.to_csv(f'{data_path}/{species}/expression_data.csv',header=False,index=False)


#--------------------------------------------------------------------------------
# # 从 CSV 文件中读取数据
data = pd.read_csv(f'{data_path}/{species}/expression_data.csv')
# print(data)
# data = data.astype(float)
# print(data.dtypes)
# 将数据写入 TSV 文件
data.to_csv(f'{data_path}/{species}/expression_data.tsv', sep='\t', header=True, index=False)


df = pd.read_csv(f'{data_path}/{species}/all_labels.csv')

# 删除 Label 为 0 的行
df = df[df['Label'] != 0]

# 将 Label 为 1 的行的 TF 和 Target 列数字加 1
#df.loc[df['Label'] == 1, ['TF', 'Target']] += 1
# 指定每列宽度为10个字符
column_widths = [10] * len(df.columns)

# 保存至txt文件，不保存列名
df.to_csv(f'{data_path}/{species}/edgelist_biogrid.txt', sep='\t',header=False, index=False)


from sklearn.preprocessing import StandardScaler

# 读取数据文件
data = pd.read_csv(f"{data_path}/{species}/expression_data.tsv", sep="\t")

# 转置数据
transposed_data = data.transpose()
print(transposed_data)
# 标准化数据
scaler = StandardScaler()
normalized_data = pd.DataFrame(scaler.fit_transform(transposed_data), columns=transposed_data.columns)

# 保存标准化后的数据到文件
normalized_data.to_csv(f"{data_path}/{species}/data_standard.txt", sep="\t", header=False)








