import pandas as pd
import pickle
from icecream import install,ic
install()

# for speci in ['sorghum']: #'wheat',,'homo','homo_C3','homo_C4''indica_BGI','japonica','zeamays''homo','indica',
    
#     datapath = f"/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/{speci}"
#     # 加载CSV文件
#     with open(f'{datapath}/label_filter_expression.pkl','rb') as fp:
#         df = pickle.load(fp)

#     #df = pd.read_csv('/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/label_sample_expression.csv', header=None)
#     ic(df)
#     # 获取第二行的数据，用于判断列的分类
#     second_row = df.iloc[0]
#     #ic(second_row)
#     # 初始化存储每个分类列索引的字典，包括第一列
#     columns_dict = {'CD': ['FEATURE_ID'], 'CK': ['FEATURE_ID'], 'CE': ['FEATURE_ID']}  # 假设第一列的索引是0

#     # 遍历第二行，根据值分类列索引，从第二列开始遍历
#     for col_index, value in second_row.items():
#         if col_index == 0:
#             continue  # 跳过第一列
#         if value in columns_dict:
#             columns_dict[value].append(col_index)
#     #ic(columns_dict)
#     # 对每个分类，提取对应的列并保存到新的CSV文件中
#     for key, col_indexes in columns_dict.items():
#         #ic(col_indexes)
#         #col_indexes = col_indexes[1:]
#         if len(col_indexes) > 1:  # 如果除了第一列之外还有其他列属于这个分类
#             # 提取列
#             #ic(df)
#             filtered_df = df.loc[:, col_indexes]
            
#             # 保存到新的CSV文件
#             filtered_df = filtered_df.iloc[1:,:]
#             filtered_df = filtered_df.dropna(axis=1)
#             ic(filtered_df.shape)
#             filtered_df.to_csv(f'{datapath}/{key}_expression.csv', index=False)

import numpy as np
for speci in ['indica_BGI','japonica','zeamays','wheat','sorghum']: #,'homo','homo_C3','homo_C4''homo','indica',
    
    datapath = f"/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count_CK_CD"
    # 加载CSV文件
    # with open(f'{datapath}/{speci}_label_sample_count.pkl','rb') as fp:
    #     df = pickle.load(fp)
    #df.columns = df.columns[1:]
    #ic(df)
    df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count_CK_CD/{speci}_label_sample_count.csv')
    ic(df)
    gene_id = df.iloc[:, [0]]  # 选择第一列作为基因ID列
    # 筛选包含 'CK' 的列（同时保留基因ID）
    ck_df = pd.concat([gene_id, df.loc[:, df.columns.str.contains('^CK', case=False)]], axis=1)

    # 筛选包含 'CD' 的列（同时保留基因ID）
    cd_df = pd.concat([gene_id, df.loc[:, df.columns.str.contains('^CD', case=False)]], axis=1)
    ic(ck_df)
    ic(cd_df)
    # 将筛选出的 DataFrame 保存到 CSV 文件中
    ck_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count_CK_CD/{speci}_CK_columns.csv', index=False)
    cd_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count_CK_CD/{speci}_CD_columns.csv', index=False)
    #df.to_csv('')
    #df = pd.read_csv('/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/label_sample_expression.csv', header=None)
    #df.index = df.index.map(str)
    
    # # 获取第二行的数据，用于判断列的分类
    # second_row = df.iloc[0]
    
    # index_list = df.index.tolist()

    # # 修改第一个元素为 'FEATURE_ID'
    # index_list[0] = 'FEATURE_ID'
    # #index_list.insert(1, np.nan)

    # # 将修改后的列表重新设置为 DataFrame 的索引
    # df.index = index_list
    # ic(df)
    # #ic(second_row)
    # # 初始化存储每个分类列索引的字典，包括第一列
    # columns_dict = {'CD': ['FEATURE_ID'], 'CK': ['FEATURE_ID'], 'CE': ['FEATURE_ID']}  # 假设第一列的索引是0

    # # 遍历第二行，根据值分类列索引，从第二列开始遍历
    # for col_index, value in second_row.items():
    #     if col_index == 0:
    #         continue  # 跳过第一列
    #     if value in columns_dict:
    #         columns_dict[value].append(col_index)
    # #ic(columns_dict)
    # # 对每个分类，提取对应的列并保存到新的CSV文件中
    # for key, col_indexes in columns_dict.items():
    #     #ic(col_indexes)
    #     #col_indexes = col_indexes[1:]
    #     if len(col_indexes) > 1:  # 如果除了第一列之外还有其他列属于这个分类
    #         # 提取列
    #         #ic(df)
    #         filtered_df = df.loc[:, col_indexes]
            
    #         # 保存到新的CSV文件
    #         filtered_df = filtered_df.iloc[1:,:]
    #         filtered_df = filtered_df.dropna(axis=1)
    #         ic(filtered_df.shape)
    #         filtered_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count_CK_CD/{speci}_{key}_count.csv', index=False)
