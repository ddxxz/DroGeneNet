
# import pickle
# import pandas as pd

# # 读取 .pkl 文件
# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/genecount_Sorghum_bicolor.pkl', 'rb') as file:
#     data = pickle.load(file)

# # 获取列名并排序
# column_names = pd.DataFrame(sorted(data.columns), columns=['Column_Names'])

# # 将排序后的列名保存到 CSV 文件，一列
# output_csv_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Sorghum_bicolor_SRRid.csv'
# column_names.to_csv(output_csv_path, index=False)

# print(f"排序后的列名已保存到 {output_csv_path}")

# import os
# import pandas as pd

# # 定义文件夹路径
# folder_path = '/home/win/4T/GeneDL/DXZ_DL/expriments/RNASeq/data/Count_wheat_drought'

# # 创建一个空的列表来存储分隔后的名字部分
# filename_list = []

# # 遍历文件夹中的所有文件
# for filename in os.listdir(folder_path):
#     # 只处理以 .tsv 结尾的文件
#     if filename.endswith('.tsv'):
#         # 以 '_' 分隔文件名，并将结果存入列表
#         name_parts = filename.split('_')
#         filename_list.append(name_parts[0])

# # 输出文件名列表
# #print(filename_list)

# # 读取 CSV 文件，假设 CSV 文件的ID列名为 'ID'
# csv_df = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Triticum_aestivum_SRRid.csv')

# # 提取 CSV 文件中的 ID 列
# csv_id_list = csv_df['Column_Names'].tolist()

# # 求交集
# intersection = set(filename_list).intersection(csv_id_list)
# print(len(intersection))
# # 求各自的特异性（不在交集中的元素）
# tsv_specific = set(filename_list) - intersection
# csv_specific = set(csv_id_list) - intersection
# print(tsv_specific)
# print(csv_specific)
# #pd.DataFrame(sorted(tsv_specific)).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/new_wheat_srr.csv',index=False)




# import pandas as pd
# sorghum_data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/CK_CD/sorghum/sorghum_CK_CD_DEG_alloutput.csv')
# wheat_data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/CK_CD/wheat/wheat_CK_CD_DEG_alloutput.csv')

# sorghum_id= sorghum_data['FEATURE_ID']
# wheat_id= wheat_data['FEATURE_ID']

# sorghum_id.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/geneid.csv',index=False)
# wheat_id.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/geneid.csv',index=False)


# import pandas as pd
# # 读取 CSV 文件
# csv_df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.txt')
# # 获取 unique 行（去重）
# unique_df = csv_df.drop_duplicates()
# print(unique_df)
# # 合并最后两列为一列
# # 假设最后两列是 'column_n-1' 和 'column_n'
# # 根据列名或者索引位置指定最后两列
# unique_df = unique_df.dropna(subset=[unique_df.columns[-2], unique_df.columns[-1]], how='all')
# merged_column = unique_df.iloc[:, -2].combine_first(unique_df.iloc[:, -1])
# # 添加到原始数据中作为一列

# unique_df['merged_column'] = merged_column

# # 输出只保留 merged_column 的结果
# result_df = unique_df[['merged_column']]
# # 查看结果
# result_df['merged_column'] = result_df['merged_column'].apply(lambda x: x.split('.')[1] if isinstance(x, str) and '.' in x else x)
# print(result_df)
# # 如果需要保存为新的 CSV 文件
# result_df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/interid.csv', index=False)

# import pandas as pd

# # 读取 Excel 文件的第一个表
# excel_file = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/Wheat_ID_Conversion.xlsx'
# df = pd.read_excel(excel_file, sheet_name=0)

# # 保留第一列和最后一列
# df = df.iloc[:, [0, -1]]

# # 删除含有 NaN 的行
# df = df.dropna()

# # 处理第一列和最后一列，保留 . 分隔的第一个值
# df.iloc[:, 0] = df.iloc[:, 0].apply(lambda x: x.split('.')[0] if isinstance(x, str) else x)
# df.iloc[:, -1] = df.iloc[:, -1].apply(lambda x: x.split('.')[0] if isinstance(x, str) else x)

# # 检查第一列是否有重复值，保留最先出现的那行
# df = df.drop_duplicates(subset=df.columns[0], keep='first')

# # 查看结果
# print(df)

# # 保存结果到新的 Excel 文件
# df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_geneid_trans.csv', index=False)

# import pandas as pd

# # 读取 txt 文件
# txt_file = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.txt'
# df_txt = pd.read_csv(txt_file)  # 假设 txt 文件是以 tab 分隔的，如果不是，请调整分隔符
# df_txt = df_txt.drop_duplicates()
# # 保留最后三列
# df_txt = df_txt.iloc[:, -3:]

# # 对倒数第二列进行处理，按照 `.` 分隔并提取第二个值作为新的列
# df_txt['new_column'] = df_txt.iloc[:, -2].apply(lambda x: x.split('.')[1] if isinstance(x, str) and '.' in x else '')

# # 读取 CSV 文件（假设该文件用于映射 Version1-Version2.1 和 UniProt-Version2.1）
# csv_file = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorgsd.csv'
# df_csv = pd.read_csv(csv_file)
# df_csv = df_csv.drop_duplicates()

# # 第一步映射：根据新的列映射到 CSV 文件中的 Version1-Version2.1
# df_txt['mapped_value1'] = df_txt['new_column'].map(df_csv.set_index('Version1')['Version2.1'])
# df_csv = df_csv.drop_duplicates(subset='UniProt', keep='first')
# df_txt['mapped_value2'] = df_txt['UniProtKB-Gene Ontology Annotation ID'].map(df_csv.set_index('UniProt')['Version2.1'])
# df_txt['merged_value'] = df_txt['mapped_value1'].combine_first(df_txt['mapped_value2'])
# # 查看最终结果
# print(df_txt)

# # 保存结果到 CSV
# df_txt.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv', index=False)

import pandas as pd
# speci='wheat'
# treat='CK'
# for layer in ['TF','PPI','KEGG']:
#     for speci in ['wheat','sorghum','homo','homo_C3','homo_C4']:
#         for treat in ['CK','CD','CE']:
#             df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{layer}/{speci}/{treat}_expression.csv')
#             print(speci,treat,df.shape)
df = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/TF_indica_BGI/CK_all.txt',sep='\t')
print(df['target_family'])

# for layer in ['TF','PPI','KEGG']:
#     for speci in ['sorghum']:
#         for treat in ['CE']:
#             file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{layer}/{speci}/{treat}_expression.csv'
#             df = pd.read_csv(file_path)

#             # 打印数据集信息
#             print(speci, treat, df.shape)

#             # 只保留前6列
#             df_subset = df.iloc[:, :6]

#             # 保存回原来的文件
#             df_subset.to_csv(file_path, index=False)


# # 读取CSV文件
# df = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/japonica/index_mapped_TF_Target.csv')

# # 假设前两列是你要分析的列
# # 生成一列来表示 A-B 和 B-A 对称的组合
# df['pair'] = df.apply(lambda row: tuple(sorted([row[0], row[1]])), axis=1)

# # 找到重复的对称组合
# duplicates = df[df.duplicated('pair', keep=False)]

# # 打印结果
# if not duplicates.empty:
#     print("以下是出现 [A, B] 和 [B, A] 对称行的数据：")
#     print(duplicates)
# else:
#     print("没有找到对称行。")
