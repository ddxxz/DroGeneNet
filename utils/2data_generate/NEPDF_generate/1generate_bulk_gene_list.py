import pandas as pd
import sys

data_type = sys.argv[1]
speci = sys.argv[2]

data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/ALL_Genefile.csv',index_col=0)
print(data)

# 对调列顺序
df = data[[data.columns[1], data.columns[0]]]  # 将第二列和第一列对调顺序

# 去除列名
df.columns = [f'Column_{i}' for i in range(len(df.columns))]  # 为列重新命名

# 保存为TXT文件
df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/NEPDF/bulk_gene_list.txt', sep='\t', index=False,header=False)









