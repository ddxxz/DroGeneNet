# import pickle
# from icecream import install,ic
# install()
# import pandas as pd
# import numpy as np
# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/zeamays_expression_batchcorrection.pkl', 'rb') as file:
#     data = pickle.load(file)
# data.reset_index(inplace=True)
# zeamay_geneid_trans=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')

# data['FEATURE_ID'] = data['FEATURE_ID'].astype(str)
# zeamay_geneid_trans['Gene stable ID'] = zeamay_geneid_trans['Gene stable ID'].astype(str)
# merged_data = pd.merge(data, zeamay_geneid_trans, how='left', left_on='FEATURE_ID', right_on='Gene stable ID')
# merged_data.drop('FEATURE_ID', axis=1, inplace=True)
# merged_data.drop('Gene stable ID', axis=1, inplace=True)
# merged_data['FEATURE_ID'] = merged_data['Oryza sativa Japonica Group gene stable ID']
# merged_data.drop('Oryza sativa Japonica Group gene stable ID', axis=1, inplace=True)
# zeamays_data = merged_data
# ic(merged_data.shape)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/wheat_expression_batchcorrection.pkl', 'rb') as file:
#     Expressiondata = pickle.load(file)
# Expressiondata.reset_index(inplace=True)
# map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
# gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
# Expressiondata['mapped_index'] = Expressiondata['FEATURE_ID'].map(gene_map)
# Expressiondata_cleaned = Expressiondata.dropna(subset=['mapped_index'])
# Expressiondata_cleaned['FEATURE_ID'] = Expressiondata_cleaned['mapped_index']
# Expressiondata_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')
# wheat_data = Expressiondata_cleaned

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/sorghum_expression_batchcorrection.pkl', 'rb') as file:
#     expression_data_df = pickle.load(file)
# expression_data_df.reset_index(inplace=True)
# map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
# gene_map = dict(zip(map_txt['Protein stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
# expression_data_df['mapped_index'] = expression_data_df['FEATURE_ID'].map(gene_map)
# expression_data_df_cleaned = expression_data_df.dropna(subset=['mapped_index'])
# expression_data_df_cleaned['FEATURE_ID'] = expression_data_df_cleaned['mapped_index']
# expression_data_df_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')
# sorghum_data = expression_data_df_cleaned

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/indica_expression_batchcorrection.pkl', 'rb') as file:
#     indica_data = pd.DataFrame(pickle.load(file))
# indica_data.reset_index(inplace=True)
# ic(indica_data.shape)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/japonica_expression_batchcorrection.pkl', 'rb') as file:
#     japonica_data = pd.DataFrame(pickle.load(file))
# ic(japonica_data.shape)
# japonica_data.reset_index(inplace=True)

# wheat_data = wheat_data.drop(columns=['SRR7185077'])
# df1=zeamays_data#.dropna()
# df2=indica_data#.dropna()
# df3=japonica_data#.dropna()
# df4=wheat_data#.dropna()
# df5=sorghum_data#.dropna()
# ic(df1,df2,df3,df4,df5)

# #common_ids = set(df1['FEATURE_ID']).intersection(df2['FEATURE_ID'], df3['FEATURE_ID'],df4['FEATURE_ID'],df5['FEATURE_ID'])
# #common_ids = set(df2['FEATURE_ID']).intersection(df3['FEATURE_ID'], df4['FEATURE_ID'])
# common_ids = set(df1['FEATURE_ID']).intersection(df5['FEATURE_ID'])
# ic(len(common_ids))
# # 从每个DataFrame中筛选包含这些共有FEATURE_ID的行
# df1_filtered = df1[df1['FEATURE_ID'].isin(common_ids)]
# df2_filtered = df2[df2['FEATURE_ID'].isin(common_ids)]
# df3_filtered = df3[df3['FEATURE_ID'].isin(common_ids)]
# df4_filtered = df4[df4['FEATURE_ID'].isin(common_ids)]
# df5_filtered = df5[df5['FEATURE_ID'].isin(common_ids)]

# df1_filtered = df1_filtered.drop_duplicates(subset='FEATURE_ID')
# df2_filtered = df2_filtered.drop_duplicates(subset='FEATURE_ID')
# df3_filtered = df3_filtered.drop_duplicates(subset='FEATURE_ID')
# df4_filtered = df4_filtered.drop_duplicates(subset='FEATURE_ID')
# df5_filtered = df5_filtered.drop_duplicates(subset='FEATURE_ID')

# #FEATURE_ID = df1_filtered['FEATURE_ID']
# FEATURE_ID = df2_filtered['FEATURE_ID']
# df2_filtered = df2_filtered.drop(columns=['FEATURE_ID'])
# df3_filtered = df3_filtered.drop(columns=['FEATURE_ID'])
# df4_filtered = df4_filtered.drop(columns=['FEATURE_ID'])
# df5_filtered = df5_filtered.drop(columns=['FEATURE_ID'])

# #merged_df = pd.concat([df1_filtered.reset_index(drop=True), df2_filtered.reset_index(drop=True), df3_filtered.reset_index(drop=True), df4_filtered.reset_index(drop=True), df5_filtered.reset_index(drop=True)], axis=1)
# merged_df = pd.concat([df2_filtered.reset_index(drop=True), df3_filtered.reset_index(drop=True), df4_filtered.reset_index(drop=True)], axis=1)
# #merged_df = pd.concat([df1_filtered.reset_index(drop=True), df5_filtered.reset_index(drop=True)], axis=1)

# #merged_df = merged_df.drop(columns=['FEATURE_ID'])
# merged_df.insert(0, 'FEATURE_ID', FEATURE_ID)

# ic(merged_df)
# nan_columns = merged_df.columns[merged_df.isna().any()]
# merged_df = merged_df.drop(columns=nan_columns[1:])
# merged_df = merged_df.dropna()

# ic(merged_df)
# merged_df.to_pickle('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/homo_C3/ExpressionData_unique_networkfilter.pkl')


#======================================================count===========================================================
import pickle
from icecream import install,ic
install()
import pandas as pd
import numpy as np
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/zeamays_count_batchcorrection.pkl', 'rb') as file:
    data = pickle.load(file)
data.reset_index(inplace=True)
zeamay_geneid_trans=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')

data['FEATURE_ID'] = data['FEATURE_ID'].astype(str)
zeamay_geneid_trans['Gene stable ID'] = zeamay_geneid_trans['Gene stable ID'].astype(str)
merged_data = pd.merge(data, zeamay_geneid_trans, how='left', left_on='FEATURE_ID', right_on='Gene stable ID')
merged_data.drop('FEATURE_ID', axis=1, inplace=True)
merged_data.drop('Gene stable ID', axis=1, inplace=True)
merged_data['FEATURE_ID'] = merged_data['Oryza sativa Japonica Group gene stable ID']
merged_data.drop('Oryza sativa Japonica Group gene stable ID', axis=1, inplace=True)
zeamays_data = merged_data
ic(merged_data.shape)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/wheat_count_batchcorrection.pkl', 'rb') as file:
    Expressiondata = pickle.load(file)
Expressiondata.reset_index(inplace=True)
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
Expressiondata['mapped_index'] = Expressiondata['FEATURE_ID'].map(gene_map)
Expressiondata_cleaned = Expressiondata.dropna(subset=['mapped_index'])
Expressiondata_cleaned['FEATURE_ID'] = Expressiondata_cleaned['mapped_index']
Expressiondata_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')
wheat_data = Expressiondata_cleaned

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/sorghum_count_batchcorrection.pkl', 'rb') as file:
    expression_data_df = pickle.load(file)
expression_data_df.reset_index(inplace=True)
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
gene_map = dict(zip(map_txt['Protein stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
expression_data_df['mapped_index'] = expression_data_df['FEATURE_ID'].map(gene_map)
expression_data_df_cleaned = expression_data_df.dropna(subset=['mapped_index'])
expression_data_df_cleaned['FEATURE_ID'] = expression_data_df_cleaned['mapped_index']
expression_data_df_cleaned.drop(columns=['mapped_index'], inplace=True, errors='ignore')
sorghum_data = expression_data_df_cleaned

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/indica_count_batchcorrection.pkl', 'rb') as file:
    indica_data = pd.DataFrame(pickle.load(file))
indica_data.reset_index(inplace=True)
ic(indica_data.shape)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/japonica_count_batchcorrection.pkl', 'rb') as file:
    japonica_data = pd.DataFrame(pickle.load(file))
ic(japonica_data.shape)
japonica_data.reset_index(inplace=True)

wheat_data = wheat_data.drop(columns=['SRR7185077'])
df1=zeamays_data#.dropna()
df2=indica_data#.dropna()
df3=japonica_data#.dropna()
df4=wheat_data#.dropna()
df5=sorghum_data#.dropna()
ic(df1,df2,df3,df4,df5)

#common_ids = set(df1['FEATURE_ID']).intersection(df2['FEATURE_ID'], df3['FEATURE_ID'],df4['FEATURE_ID'],df5['FEATURE_ID'])
#common_ids = set(df2['FEATURE_ID']).intersection(df3['FEATURE_ID'], df4['FEATURE_ID'])
common_ids = set(df1['FEATURE_ID']).intersection(df5['FEATURE_ID'])
ic(len(common_ids))
# 从每个DataFrame中筛选包含这些共有FEATURE_ID的行
df1_filtered = df1[df1['FEATURE_ID'].isin(common_ids)]
df2_filtered = df2[df2['FEATURE_ID'].isin(common_ids)]
df3_filtered = df3[df3['FEATURE_ID'].isin(common_ids)]
df4_filtered = df4[df4['FEATURE_ID'].isin(common_ids)]
df5_filtered = df5[df5['FEATURE_ID'].isin(common_ids)]

df1_filtered = df1_filtered.drop_duplicates(subset='FEATURE_ID')
df2_filtered = df2_filtered.drop_duplicates(subset='FEATURE_ID')
df3_filtered = df3_filtered.drop_duplicates(subset='FEATURE_ID')
df4_filtered = df4_filtered.drop_duplicates(subset='FEATURE_ID')
df5_filtered = df5_filtered.drop_duplicates(subset='FEATURE_ID')

#FEATURE_ID = df1_filtered['FEATURE_ID']
FEATURE_ID = df2_filtered['FEATURE_ID']
df2_filtered = df2_filtered.drop(columns=['FEATURE_ID'])
df3_filtered = df3_filtered.drop(columns=['FEATURE_ID'])
df4_filtered = df4_filtered.drop(columns=['FEATURE_ID'])
df5_filtered = df5_filtered.drop(columns=['FEATURE_ID'])

#merged_df = pd.concat([df1_filtered.reset_index(drop=True), df2_filtered.reset_index(drop=True), df3_filtered.reset_index(drop=True), df4_filtered.reset_index(drop=True), df5_filtered.reset_index(drop=True)], axis=1)
merged_df = pd.concat([df2_filtered.reset_index(drop=True), df3_filtered.reset_index(drop=True), df4_filtered.reset_index(drop=True)], axis=1)
#merged_df = pd.concat([df1_filtered.reset_index(drop=True), df5_filtered.reset_index(drop=True)], axis=1)

#merged_df = merged_df.drop(columns=['FEATURE_ID'])
merged_df.insert(0, 'FEATURE_ID', FEATURE_ID)

ic(merged_df)
nan_columns = merged_df.columns[merged_df.isna().any()]
merged_df = merged_df.drop(columns=nan_columns[1:])
merged_df = merged_df.dropna()

ic(merged_df)
merged_df.to_pickle('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/homo_C3_count_batchcorrection.pkl')
