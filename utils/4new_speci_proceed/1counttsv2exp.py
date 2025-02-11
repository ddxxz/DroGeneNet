#-------------------------------------拼接各个tsv脚本---------------------------------------------
# import pandas as pd
# import os
# from icecream import install,ic
# install()
# from tqdm import tqdm
# # 指定你的文件夹路径
# folder_path = '/home/win/4T/GeneDL/DXZ_DL/expriments/RNASeq/data/Count_sorghum_drought'
# output_dir ='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci'
# species = 'Sorghum_bicolor'
# csv_file = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Sorghum_bicolor_SRRid.csv'
# csv_df = pd.read_csv(csv_file)
# srr_id_list = csv_df['Column_Names'].tolist()  # 提取CSV中的ID列

# # 初始化用于合并所有文件的DataFrame
# merged_df = pd.DataFrame()

# # 遍历指定文件夹查找TSV文件
# for filename in tqdm(os.listdir(folder_path)):
#     if filename.endswith('.tsv'):
#         # 提取文件名中的SRRid部分，去掉'_featureCounts'等后缀
#         srr_id = filename.split('.')[0].replace('_featureCounts', '')
        
#         # 如果SRRid不在CSV文件中的ID列表中，则跳过该文件
#         if srr_id not in srr_id_list:
#             continue
        
#         # 构建文件路径
#         file_path = os.path.join(folder_path, filename)
        
#         # 读取TSV文件，并设置列名
#         df = pd.read_csv(file_path, sep='\t', skiprows=1, header=None)
#         df.columns = ['Geneid', 'Unknown1', 'Start', 'End', 'Strand', 'Length', 'Count']
        
#         # 文件名命名新列
#         new_column_name = srr_id
        
#         # 如果是第一个文件，初始化merged_df
#         if merged_df.empty:
#             merged_df['Geneid'] = df['Geneid']
#             merged_df['Length'] = df['Length']
#             merged_df[new_column_name] = df['Count']  # 添加第七列（Count列）
#         else:
#             # 将后续文件的Count列添加到merged_df中，匹配Geneid
#             merged_df[new_column_name] = df['Count']
# # 检查结果
# print(merged_df)
# merged_df  = merged_df.iloc[1:,:]
# output_path = f'{output_dir}/genecount_{species}.pkl'
# # 使用to_pickle方法保存DataFrames
# merged_df.to_pickle(output_path)

# #-------------------------------------生成表达矩阵和基因长度文件------------------------------------
# import pandas as pd

# # 读取CSV文件
# df = pd.read_pickle(output_path) #pd.read_csv('/mnt/h/study/deep_learning/gene/RNA_Seq/genecount/drought_samplecount.csv')

# # 创建一个新的DataFrame来保存gene_id和gene_length
# lengths_df = df.iloc[:, :2]  # 选择前两列
# lengths_df.columns = ['FEATURE_ID', 'GENE_LENGTHS']  # 更改列名

# # 保存第一列和第二列为新的CSV文件
# lengths_df.to_pickle(f'{output_dir}/lengths.pkl')

# # 为了转置剩余的列（除了第二列外），先删除第二列
# exp_df = df.drop(df.columns[1], axis=1)  # 删除第二列

# # 转置剩余列
# #exp_df_transposed = exp_df.T  # 只转置第二列之后的数据
# exp_df.rename(columns={exp_df.columns[0]: 'FEATURE_ID'}, inplace=True)
# # 保存转置后的DataFrame到CSV
# print(exp_df)
# feature_id = exp_df.iloc[:,0]
# exp_df = exp_df.iloc[:,1:]

# exp_df.to_pickle(f'{output_dir}/exp.pkl')
# #ic(exp_df)
# #-------------------------------------标准化------------------------------------

# import pandas as pd
# from icecream import install,ic
# install()
# # Load the expression counts data (assuming it's in CSV format)
# # The DataFrame should have genes as rows and samples as columns
# exp_counts_df = pd.read_pickle(f'{output_dir}/exp.pkl')
# ic(exp_counts_df)
# #exp_counts_df = exp_counts_df.iloc[:,1:]
# # Load the gene lengths data
# lengths_df = pd.read_pickle(f'{output_dir}/lengths.pkl')
# lengths_df = lengths_df.iloc[:,1:]
# #ic(lengths_df)
# #ic(exp_counts_df.shape,lengths_df.shape,lengths_df.values.squeeze().shape)
# # --- 计算 CPM ---
# # Sum the counts for each sample (column) to get the total counts
# exp_counts_df = exp_counts_df.apply(pd.to_numeric, errors='coerce')

# # 检查 exp_counts_df 中是否有 NaN 值
# print("Checking for NaN values in exp_counts_df:")
# print(exp_counts_df.isna().sum())
# # 填充 NaN 值为 0
# exp_counts_df = exp_counts_df.fillna(0)
# # 检查数据框的形状，确保数据存在
# print("Shape of exp_counts_df:", exp_counts_df.shape)

# total_counts_per_sample = exp_counts_df.sum(axis=0)
# ic(total_counts_per_sample)

# total_counts_per_sample = total_counts_per_sample.replace(0, 1e-6).fillna(1e-6)
# ic(total_counts_per_sample)
# # Divide each count by the total counts and multiply by one million
# cpm_df = exp_counts_df.divide(total_counts_per_sample, axis=1) * 10**6
# ic(cpm_df)
# # --- 计算 FPKM ---
# # Divide the counts by the lengths of the genes to get counts per kilobase
# # lengths_df.values.squeeze() converts the DataFrame to a single column array
# # The lengths are assumed to be in base pairs, so we divide by 1,000 to convert to kilobases
# ic(exp_counts_df)
# ic(lengths_df)
# # 显示数据类型
# ic(lengths_df.dtypes)
# # 输出所有可能非数值的行
# ic(lengths_df[pd.to_numeric(lengths_df.squeeze(), errors='coerce').isna()])
# lengths_df = lengths_df.apply(pd.to_numeric, errors='coerce')
# # 再次检查数据
# ic(lengths_df.dtypes)
# exp_counts_df = exp_counts_df.apply(pd.to_numeric, errors='coerce')
# counts_per_kb = exp_counts_df.divide(lengths_df['GENE_LENGTHS'] / 10**3, axis=0)
# ic(counts_per_kb)
# # Then divide by the total counts (in millions) to get FPKM
# fpkm_df = counts_per_kb.divide(total_counts_per_sample / 10**6 +1e-6, axis=1)
# ic(fpkm_df)
# # --- 计算 TPM ---
# # First, get the per million scaling factor
# # Divide counts per kb (from the FPKM calculation step) by the total counts per sample
# per_million_scaling_factor = counts_per_kb.divide(total_counts_per_sample, axis=1).sum(axis=0)
# ic(per_million_scaling_factor)
# # Divide counts per kb by the per million scaling factor, and then normalize by multiplying by one million
# tpm_df = counts_per_kb.divide(per_million_scaling_factor, axis=1) * 10**6
# ic(tpm_df)
# # Save the CPM, FPKM, and TPM data to new CSV files
# cpm_df_with_id = pd.concat([feature_id, cpm_df], axis=1)
# fpkm_df_with_id = pd.concat([feature_id, fpkm_df], axis=1)
# tpm_df_with_id = pd.concat([feature_id, tpm_df], axis=1)
# cpm_df_with_id.to_pickle(f'{output_dir}/exp_cpm_transcript_{species}.pkl')
# fpkm_df_with_id.to_pickle(f'{output_dir}/exp_fpkm_transcript_{species}.pkl')
# tpm_df_with_id.to_pickle(f'{output_dir}/exp_tpm_transcript_{species}.pkl')

#-------------------------------------生成表达矩阵和基因长度文件------------------------------------


# import pickle

# # 读取 .pkl 文件
# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/genecount_Triticum_aestivum.pkl', 'rb') as file:
#     data = pickle.load(file)

# data = data.drop(columns=['Length'])  # 删除 Length 列
# data = data.rename(columns={'Geneid': 'FEATURE_ID'})  # 将 Geneid 列重命名为 FEATURE_ID

# # 查看修改后的数据（可选）
# print(data.head())

# # 保存修改后的数据回到 .pkl 文件
# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/genecount_Triticum_aestivum.pkl', 'wb') as file:
#     pickle.dump(data, file)
# # 查看读取的数据
# print(data)


import pickle
import pandas as pd
from tqdm import tqdm
# 读取 .pkl 文件
species = 'Triticum_aestivum'
with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/exp_cpm_transcript_Triticum_aestivum.pkl', 'rb') as file:
    data = pickle.load(file)
# 提取Geneid的基准列，按`-`符号分割，并取第一个值
merged_df=data
# merged_df['Geneid_base'] = merged_df['FEATURE_ID'].apply(lambda x: x.split('-')[0])
print(merged_df)

def extract_base_id(feature_id):
    # 先去除末尾的 '-E1' 或类似的后缀
    if '-' in feature_id and feature_id.endswith('-E1'):
        feature_id = feature_id.split('-')[0]
    
    # 然后处理 '.' 分隔的部分
    if '.' in feature_id:
        return feature_id.split('.')[0]
    else:
        return feature_id

# 应用函数提取基因基础部分
merged_df['Geneid_base'] = merged_df['FEATURE_ID'].apply(extract_base_id)
print(merged_df)

# 确保所有需要求和的列是数值类型
# for col in tqdm(merged_df.columns):
#     if col not in ['FEATURE_ID', 'Geneid_base']:
#         merged_df[col] = pd.to_numeric(merged_df[col], errors='coerce')  # 转换为数值类型
# 直接对非字符列进行批量转换
from joblib import Parallel, delayed
import pandas as pd

def process_column(col):
    return pd.to_numeric(col, errors='coerce')

# 获取需要转换的数值列
numeric_columns = merged_df.columns.difference(['FEATURE_ID', 'Geneid_base'])
# 并行处理每列，并将结果转换为 DataFrame
processed_columns = Parallel(n_jobs=-1)(delayed(process_column)(merged_df[col]) for col in numeric_columns)
# 将并行处理的结果转换回 DataFrame
processed_df = pd.DataFrame(processed_columns).T  # 需要转置，因为并行化处理返回的是按列的数据
# 将转换后的结果重新赋值回原来的 DataFrame
merged_df[numeric_columns] = processed_df.values
print(merged_df)

# 根据Geneid_base进行分组，并对其数值列（除了Geneid和Length之外）进行汇总
grouped_df = merged_df.groupby('Geneid_base').agg({
    #'Length': 'first',  # Length列保留第一条（或可选择汇总方式）
    **{col: 'sum' for col in merged_df.columns if col not in ['FEATURE_ID', 'Geneid_base']}
}).reset_index()

# 删除辅助列 Geneid_base
#grouped_df.drop(columns=['Geneid_base'], inplace=True)

# 检查结果
print(grouped_df)

# 去掉第一行示例代码行
grouped_df = grouped_df.iloc[1:, :]
grouped_df = grouped_df.rename(columns={'Geneid_base': 'FEATURE_ID'})
# 使用to_pickle方法保存DataFrames
output_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/gene_cpm_{species}.pkl'
grouped_df.to_pickle(output_path)
