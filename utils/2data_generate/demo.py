
import pandas as pd

# 读取CSV文件
file_path = '/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/japonica/CK_expression.csv'  # 替换为实际的文件路径
data = pd.read_csv(file_path)

nan_counts = data.isnull().sum()

# 打印包含NaN值的列名和对应的NaN值数量
nan_columns = nan_counts[nan_counts > 0]

if not nan_columns.empty:
    print("The CSV file contains NaN values in the following columns:")
    print(nan_columns)
else:
    print("The CSV file does not contain any NaN values.")

print(data['SRR12764688'])










