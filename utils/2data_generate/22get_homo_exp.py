import pickle
from icecream import install,ic
install()
import pandas as pd
import numpy as np
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/zeamays_expression_batchcorrection.pkl', 'rb') as file:
    data = pickle.load(file)
data.reset_index(inplace=True)
zeamay_geneid_trans=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')

#data['FEATURE_ID'].to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid.csv',index=False)

data['FEATURE_ID'] = data['FEATURE_ID'].astype(str)
zeamay_geneid_trans['Gene stable ID'] = zeamay_geneid_trans['Gene stable ID'].astype(str)

merged_data = pd.merge(data, zeamay_geneid_trans, how='left', left_on='FEATURE_ID', right_on='Gene stable ID')
merged_data.drop('FEATURE_ID', axis=1, inplace=True)
merged_data.drop('Gene stable ID', axis=1, inplace=True)
merged_data['FEATURE_ID'] = merged_data['Oryza sativa Japonica Group gene stable ID']
merged_data.drop('Oryza sativa Japonica Group gene stable ID', axis=1, inplace=True)

ic(merged_data.shape)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/indica_expression_batchcorrection.pkl', 'rb') as file:
    indica_data = pd.DataFrame(pickle.load(file))
indica_data.reset_index(inplace=True)
ic(indica_data.shape)
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/japonica_expression_batchcorrection.pkl', 'rb') as file:
    japonica_data = pd.DataFrame(pickle.load(file))
ic(japonica_data.shape)
japonica_data.reset_index(inplace=True)
df1=merged_data
df2=indica_data
df3=japonica_data


common_ids = set(df1['FEATURE_ID']).intersection(df2['FEATURE_ID'], df3['FEATURE_ID'])
ic(len(common_ids))
# 从每个DataFrame中筛选包含这些共有FEATURE_ID的行
df1_filtered = df1[df1['FEATURE_ID'].isin(common_ids)]
df2_filtered = df2[df2['FEATURE_ID'].isin(common_ids)]
df3_filtered = df3[df3['FEATURE_ID'].isin(common_ids)]

df1_filtered = df1_filtered.drop_duplicates(subset='FEATURE_ID')
df2_filtered = df2_filtered.drop_duplicates(subset='FEATURE_ID')
df3_filtered = df3_filtered.drop_duplicates(subset='FEATURE_ID')

def move_last_column_to_first(df):
    # 获取列名列表
    cols = df.columns.tolist()
    # 将最后一列名移到列表的开头
    cols = [cols[-1]] + cols[:-1]
    # 重新排序DataFrame的列
    new_df = df[cols]
    return new_df

# 对 df1_filtered, df2_filtered, df3_filtered 应用这个函数
df1_filtered = move_last_column_to_first(df1_filtered)
# df2_filtered = move_last_column_to_first(df2_filtered)
# df3_filtered = move_last_column_to_first(df3_filtered)

#ic(df1_filtered.shape,df2_filtered.shape,df3_filtered.shape)
df1_filtered = df1_filtered.set_index('FEATURE_ID').reindex(df3_filtered['FEATURE_ID']).reset_index()
df2_filtered = df2_filtered.set_index('FEATURE_ID').reindex(df3_filtered['FEATURE_ID']).reset_index()
#ic(df1_filtered,df2_filtered,df3_filtered)

final_array = np.concatenate([df1_filtered.values, df2_filtered.drop(columns=['FEATURE_ID']).values, df3_filtered.drop(columns=['FEATURE_ID']).values], axis=1)

# 构建新的 DataFrame，保留列名
final_columns = ['FEATURE_ID'] + [col for col in df1_filtered.columns if col != 'FEATURE_ID'] + list(df2_filtered.columns[1:]) + list(df3_filtered.columns[1:])
final_df = pd.DataFrame(final_array, columns=final_columns)

# 将 FEATURE_ID 列移动到第一列
final_df = final_df[['FEATURE_ID'] + [col for col in final_df.columns if col != 'FEATURE_ID']]

# 拼接三个表格，并只保留一列 FEATURE_ID
#final_df = pd.concat([df3_filtered['FEATURE_ID'], df1_filtered.drop(columns=['FEATURE_ID']), df2_filtered.drop(columns=['FEATURE_ID']), df3_filtered.drop(columns=['FEATURE_ID'])], axis=1,ignore_index=True)
# 将 FEATURE_ID 列移动到第一列
#final_df = final_df[['FEATURE_ID'] + [col for col in final_df.columns if col != 'FEATURE_ID']]
# # 拼接这些DataFrame
# ic(len(df1_filtered['FEATURE_ID'].unique()))
# ic(df2_filtered)
# ic(df3_filtered)
# final_df = pd.concat([df1_filtered, df2_filtered, df3_filtered], axis=1,ignore_index=False)

filtered_exp_df= final_df.loc[:, ~final_df.columns.duplicated()]
#ic(final_df)


# zeamays_tf_target_data=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/zeamays/TF_Target.csv')

# mapping_dict = dict(zip(zeamay_geneid_trans['Gene stable ID'], zeamay_geneid_trans['Oryza sativa Japonica Group gene stable ID']))

# # 使用映射字典对 df2 的两列进行映射
# zeamays_tf_target_data['TF'] = zeamays_tf_target_data['TF'].map(mapping_dict)
# zeamays_tf_target_data['Target'] = zeamays_tf_target_data['Target'].map(mapping_dict)
# #ic(zeamays_tf_target_data)


# indica_tf_target_data=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/indica/TF_Target.csv')
# japonica_tf_target_data=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/japonica/TF_Target.csv')

# merged_df = indica_tf_target_data.merge(japonica_tf_target_data, on=['TF', 'Target']).merge(zeamays_tf_target_data, on=['TF', 'Target'])
# #ic(merged_df)
# gene=pd.concat([merged_df['TF'],merged_df['Target']],axis=0)
# gene_unique = gene.unique()

# ic(len(gene.unique()))

# intersection = list(set(final_df['FEATURE_ID'].tolist()) & set(gene_unique.tolist()))
# ic(len(intersection))

# filtered_exp_df = final_df[final_df['FEATURE_ID'].isin(gene_unique)]
# final_tf_target=merged_df
# ic(filtered_exp_df)
# ic(final_tf_target)

filtered_exp_df.to_pickle('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/homo/ExpressionData_unique_networkfilter.pkl')
#final_tf_target.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/TF_Target.csv',index=False)


#=====================================================genecount==========================================================================
# import pickle
# from icecream import install,ic
# install()
# import pandas as pd
# import numpy as np
# with open('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/genecount_gene_zeamays.pkl', 'rb') as file:
#     data = pickle.load(file)

# zeamay_geneid_trans=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')

# #data['FEATURE_ID'].to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid.csv',index=False)

# data['FEATURE_ID'] = data['FEATURE_ID'].astype(str)
# zeamay_geneid_trans['Gene stable ID'] = zeamay_geneid_trans['Gene stable ID'].astype(str)

# merged_data = pd.merge(data, zeamay_geneid_trans, how='left', left_on='FEATURE_ID', right_on='Gene stable ID')
# merged_data.drop('FEATURE_ID', axis=1, inplace=True)
# merged_data.drop('Gene stable ID', axis=1, inplace=True)
# merged_data['FEATURE_ID'] = merged_data['Oryza sativa Japonica Group gene stable ID']
# merged_data.drop('Oryza sativa Japonica Group gene stable ID', axis=1, inplace=True)

# ic(merged_data.shape)

# with open('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/genecount_gene_OS_Oryza_indica.pkl', 'rb') as file:
#     indica_data = pd.DataFrame(pickle.load(file))
# ic(indica_data.shape)
# with open('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/japonica/genecount_gene_japonica.pkl', 'rb') as file:
#     japonica_data = pd.DataFrame(pickle.load(file))
# ic(japonica_data.shape)

# df1=merged_data
# df2=indica_data
# df3=japonica_data


# common_ids = set(df1['FEATURE_ID']).intersection(df2['FEATURE_ID'], df3['FEATURE_ID'])
# ic(len(common_ids))
# # 从每个DataFrame中筛选包含这些共有FEATURE_ID的行
# df1_filtered = df1[df1['FEATURE_ID'].isin(common_ids)]
# df2_filtered = df2[df2['FEATURE_ID'].isin(common_ids)]
# df3_filtered = df3[df3['FEATURE_ID'].isin(common_ids)]

# df1_filtered = df1_filtered.drop_duplicates(subset='FEATURE_ID')
# df2_filtered = df2_filtered.drop_duplicates(subset='FEATURE_ID')
# df3_filtered = df3_filtered.drop_duplicates(subset='FEATURE_ID')

# def move_last_column_to_first(df):
#     # 获取列名列表
#     cols = df.columns.tolist()
#     # 将最后一列名移到列表的开头
#     cols = [cols[-1]] + cols[:-1]
#     # 重新排序DataFrame的列
#     new_df = df[cols]
#     return new_df

# # 对 df1_filtered, df2_filtered, df3_filtered 应用这个函数
# df1_filtered = move_last_column_to_first(df1_filtered)
# # df2_filtered = move_last_column_to_first(df2_filtered)
# # df3_filtered = move_last_column_to_first(df3_filtered)

# #ic(df1_filtered.shape,df2_filtered.shape,df3_filtered.shape)
# df1_filtered = df1_filtered.set_index('FEATURE_ID').reindex(df3_filtered['FEATURE_ID']).reset_index()
# df2_filtered = df2_filtered.set_index('FEATURE_ID').reindex(df3_filtered['FEATURE_ID']).reset_index()
# #ic(df1_filtered,df2_filtered,df3_filtered)

# final_array = np.concatenate([df1_filtered.values, df2_filtered.drop(columns=['FEATURE_ID']).values, df3_filtered.drop(columns=['FEATURE_ID']).values], axis=1)

# # 构建新的 DataFrame，保留列名
# final_columns = ['FEATURE_ID'] + [col for col in df1_filtered.columns if col != 'FEATURE_ID'] + list(df2_filtered.columns[1:]) + list(df3_filtered.columns[1:])
# final_df = pd.DataFrame(final_array, columns=final_columns)

# # 将 FEATURE_ID 列移动到第一列
# final_df = final_df[['FEATURE_ID'] + [col for col in final_df.columns if col != 'FEATURE_ID']]

# # 拼接三个表格，并只保留一列 FEATURE_ID
# #final_df = pd.concat([df3_filtered['FEATURE_ID'], df1_filtered.drop(columns=['FEATURE_ID']), df2_filtered.drop(columns=['FEATURE_ID']), df3_filtered.drop(columns=['FEATURE_ID'])], axis=1,ignore_index=True)
# # 将 FEATURE_ID 列移动到第一列
# #final_df = final_df[['FEATURE_ID'] + [col for col in final_df.columns if col != 'FEATURE_ID']]
# # # 拼接这些DataFrame
# # ic(len(df1_filtered['FEATURE_ID'].unique()))
# # ic(df2_filtered)
# # ic(df3_filtered)
# # final_df = pd.concat([df1_filtered, df2_filtered, df3_filtered], axis=1,ignore_index=False)

# final_df= final_df.loc[:, ~final_df.columns.duplicated()]
# ic(final_df)
# final_df.to_pickle('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/homo_genecount.pkl')


