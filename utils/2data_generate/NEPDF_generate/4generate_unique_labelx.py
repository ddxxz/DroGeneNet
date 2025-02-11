import pandas as pd
import sys
data_type = sys.argv[1]
speci = sys.argv[2]
# 读取CSV文件
df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/all_labels.csv')

# 创建新的DataFrame
new_df = pd.DataFrame(columns=['TF', 'Target', 'Label'])
label1_rows = df[df['Label'] == 1]
label0_rows = df[df['Label'] == 0]

# 确保两个DataFrame的行数一致
min_len = min(len(label1_rows), len(label0_rows))
label1_rows = label1_rows.iloc[:min_len]
label0_rows = label0_rows.iloc[:min_len]

# 合并Label为1和Label为0的行
new_df = pd.concat([label1_rows, label0_rows]).sort_index(kind='merge')
# # 提取Label为1的行
# while len(df[df['Label'] == 1]) > 0:
#     # 提取Label为1的行
#     label1_row = df[df['Label'] == 1].iloc[0]

#     # # 复制Label为1的行，设置标记为2
#     # new_row = label1_row.copy()
#     # new_row['Label'] = 2
#     # new_row['TF'] = label1_row['Target']
#     # new_row['Target'] = label1_row['TF']
#     # # 将新行添加到新的DataFrame中
#     new_df = pd.concat([new_df, label1_row.to_frame().T])
#     #new_df = pd.concat([new_df, new_row.to_frame().T])

#     # 提取Label为0的行
#     label0_row = df[df['Label'] == 0].iloc[0]

#     # 将Label为0的行插入到新的DataFrame的下一个位置
#     new_df = pd.concat([new_df, label0_row.to_frame().T])

#     # 删除已处理的行
#     df = df.drop([label1_row.name, label0_row.name])

# 保存为TXT文件并不包含列名
new_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/NEPDF/unique_labelx.txt', sep='\t', index=False, header=False)

print("数据已保存到output.txt文件中。")