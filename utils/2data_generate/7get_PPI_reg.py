# import pandas as pd

# # 读取txt文件
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica.39947.protein.links.v12.0.txt'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep=' ', header=0)
# print(data)
# # 去掉protein1和protein2列中的数字部分
# data['protein1'] = data['protein1'].apply(lambda x: x.split('.')[1])
# data['protein2'] = data['protein2'].apply(lambda x: x.split('.')[1])

# # 提取需要的列
# result = data[['protein1', 'protein2','combined_score']]
# filtered_result = result[result['combined_score'] > 700]
# print(filtered_result)

# filtered_result['sorted_proteins'] = filtered_result.apply(lambda x: tuple(sorted([x['protein1'], x['protein2']])), axis=1)

# # 去除重复的行，只保留唯一的蛋白质对
# filtered_result = filtered_result.drop_duplicates(subset='sorted_proteins').drop('sorted_proteins', axis=1)
# print(filtered_result)
# # 保存处理后的数据至csv文件
# output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/PPI_reg.csv'  # 替换为实际的输出文件路径
# filtered_result.to_csv(output_file_path, index=False)

# # 合并两列并获取唯一值
# unique_proteins = pd.unique(filtered_result[['protein1', 'protein2']].values.ravel())

# # 保存唯一蛋白质ID至另一个csv文件
# unique_output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica_unique_proteins.csv'  # 替换为实际的输出文件路径
# pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)

# print("Processed data saved to:", output_file_path)
# print("Unique proteins saved to:", unique_output_file_path)



# import pandas as pd

# # 读取txt文件
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/zeamays.4577.protein.links.v12.0.txt'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep=' ', header=0)
# print(data)
# # 去掉protein1和protein2列中的数字部分
# data['protein1'] = data['protein1'].apply(lambda x: x.split('.')[1])
# data['protein2'] = data['protein2'].apply(lambda x: x.split('.')[1])

# # 提取需要的列
# result = data[['protein1', 'protein2','combined_score']]
# filtered_result = result[result['combined_score'] > 700]
# print(filtered_result)

# filtered_result['sorted_proteins'] = filtered_result.apply(lambda x: tuple(sorted([x['protein1'], x['protein2']])), axis=1)

# # 去除重复的行，只保留唯一的蛋白质对
# filtered_result = filtered_result.drop_duplicates(subset='sorted_proteins').drop('sorted_proteins', axis=1)
# print(filtered_result)
# # 保存处理后的数据至csv文件
# output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/PPI_reg.csv'  # 替换为实际的输出文件路径
# filtered_result.to_csv(output_file_path, index=False)

# # 合并两列并获取唯一值
# unique_proteins = pd.unique(filtered_result[['protein1', 'protein2']].values.ravel())

# # 保存唯一蛋白质ID至另一个csv文件
# unique_output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/zeamays_unique_proteins.csv'  # 替换为实际的输出文件路径
# pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)

# print("Processed data saved to:", output_file_path)
# print("Unique proteins saved to:", unique_output_file_path)



# import pandas as pd
# from icecream import install,ic
# install()

# # 读取txt文件
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/indica.39946.protein.links.v12.0.txt'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep=' ', header=0)
# print(data)
# # 去掉protein1和protein2列中的数字部分
# data['protein1'] = data['protein1'].apply(lambda x: x.split('.')[1])
# data['protein2'] = data['protein2'].apply(lambda x: x.split('.')[1])

# # 提取需要的列
# result = data[['protein1', 'protein2','combined_score']]
# filtered_result = result[result['combined_score'] > 700]
# print(filtered_result)

# # df_reversed = filtered_result.rename(columns={"protein1": "protein2", "protein2": "protein1"})
# # # 将原始DataFrame与修改后的DataFrame进行合并，查找两个方向都存在的对
# # mutual_regulation = pd.merge(filtered_result, df_reversed, on=["protein1", "protein2"])
# # # 查看结果
# # mutual_regulation_count = len(mutual_regulation) // 2  # 因为结果中每对会出现两次，所以需要除以2
# # ic(mutual_regulation_count)
# filtered_result['sorted_proteins'] = filtered_result.apply(lambda x: tuple(sorted([x['protein1'], x['protein2']])), axis=1)

# # 去除重复的行，只保留唯一的蛋白质对
# filtered_result = filtered_result.drop_duplicates(subset='sorted_proteins').drop('sorted_proteins', axis=1)
# print(filtered_result)

# # 保存处理后的数据至csv文件
# output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/PPI_reg.csv'  # 替换为实际的输出文件路径
# filtered_result.to_csv(output_file_path, index=False)

# # 合并两列并获取唯一值
# unique_proteins = pd.unique(filtered_result[['protein1', 'protein2']].values.ravel())

# # 保存唯一蛋白质ID至另一个csv文件
# unique_output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/indica_BGI_unique_proteins.csv'  # 替换为实际的输出文件路径
# pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)

# print("Processed data saved to:", output_file_path)
# print("Unique proteins saved to:", unique_output_file_path)



# import pandas as pd

# # 读取txt文件
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/sorghum.4558.protein.links.v12.0.txt'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep=' ', header=0)
# print(data)
# # 去掉protein1和protein2列中的数字部分
# data['protein1'] = data['protein1'].apply(lambda x: x.split('.')[1])
# data['protein2'] = data['protein2'].apply(lambda x: x.split('.')[1])

# # 提取需要的列
# result = data[['protein1', 'protein2','combined_score']]
# filtered_result = result[result['combined_score'] > 700]
# print(filtered_result)

# filtered_result['sorted_proteins'] = filtered_result.apply(lambda x: tuple(sorted([x['protein1'], x['protein2']])), axis=1)

# # 去除重复的行，只保留唯一的蛋白质对
# filtered_result = filtered_result.drop_duplicates(subset='sorted_proteins').drop('sorted_proteins', axis=1)
# print(filtered_result)
# # 保存处理后的数据至csv文件
# output_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/PPI_reg.csv'  # 替换为实际的输出文件路径
# filtered_result.to_csv(output_file_path, index=False)

# # 合并两列并获取唯一值
# unique_proteins = pd.unique(filtered_result[['protein1', 'protein2']].values.ravel())

# # 保存唯一蛋白质ID至另一个csv文件
# unique_output_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/sorghum/sorghum_unique_proteins.csv'  # 替换为实际的输出文件路径
# pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)

# print("Processed data saved to:", output_file_path)
# print("Unique proteins saved to:", unique_output_file_path)





import pandas as pd
from tqdm import tqdm

# 读取txt文件，设置 chunksize 来逐块读取文件
file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/wheat/wheat.4565.protein.links.v12.0.txt'  # 替换为实际的文件路径
chunksize = 100000  # 每次读取 100,000 行

# 获取文件的总行数
total_lines = sum(1 for line in open(file_path, 'r'))

# 计算总块数
total_chunks = total_lines // chunksize + 1

# 初始化空的DataFrame，用于存储过滤后的结果
filtered_result = pd.DataFrame()

# 逐块读取并处理数据，使用tqdm显示进度条
for chunk in tqdm(pd.read_csv(file_path, sep=' ', header=0, chunksize=chunksize), total=total_chunks):
    # 去掉protein1和protein2列中的数字部分
    chunk['protein1'] = chunk['protein1'].apply(lambda x: x.split('.')[1])
    chunk['protein2'] = chunk['protein2'].apply(lambda x: x.split('.')[1])
    
    # 提取需要的列，并过滤combined_score > 700的行
    result = chunk[['protein1', 'protein2', 'combined_score']]
    chunk_filtered_result = result[result['combined_score'] > 700]
    
    # 对 protein1 和 protein2 进行排序，创建新的列 sorted_proteins
    chunk_filtered_result['sorted_proteins'] = chunk_filtered_result.apply(
        lambda x: tuple(sorted([x['protein1'], x['protein2']])), axis=1)
    
    # 去除重复的行，只保留唯一的蛋白质对
    chunk_filtered_result = chunk_filtered_result.drop_duplicates(subset='sorted_proteins').drop('sorted_proteins', axis=1)
    
    # 将当前块的过滤结果追加到最终结果 DataFrame
    filtered_result = pd.concat([filtered_result, chunk_filtered_result])

# 保存处理后的数据至csv文件
output_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/wheat/PPI_reg.csv'  # 替换为实际的输出文件路径
filtered_result.to_csv(output_file_path, index=False)

# 合并两列并获取唯一值
unique_proteins = pd.unique(filtered_result[['protein1', 'protein2']].values.ravel())

# 保存唯一蛋白质ID至另一个csv文件
unique_output_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/wheat/wheat_unique_proteins.csv'  # 替换为实际的输出文件路径
pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)






