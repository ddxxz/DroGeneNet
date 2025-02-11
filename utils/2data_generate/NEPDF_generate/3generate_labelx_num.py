import math
import sys
data_type = sys.argv[1]
speci = sys.argv[2]

# 读取txt文件并计算行数
file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/NEPDF/bulk_gene_list.txt'  # 替换为你的txt文件路径
with open(file_path, 'r', encoding='utf-8') as file:
    lines = file.readlines()
    num_lines = len(lines)

# 计算每份的行数，最后一份可能会少一些
lines_per_group = math.ceil(num_lines / 10)

# 生成平分区间开始的索引
indices = [i * lines_per_group for i in range(10)]
# 调整最后一个索引，不超过总行数
indices[-1] = max(indices[-1], num_lines - lines_per_group)

# 将结果输出到新的txt文件中
output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/NEPDF/labelx_num.txt'  # 替换为输出文件的路径
with open(output_file_path, 'w', encoding='utf-8') as output_file:
    for index in indices:
        output_file.write(f"{index}\n")

print("分区索引已成功保存到文件。")