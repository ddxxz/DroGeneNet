
import pickle
from icecream import install,ic
install()
species='wheat' #sorghum
with open(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/gene_fpkm_Triticum_aestivum.pkl','rb') as fp:
    data=pickle.load(fp)
print(data)
data.set_index('FEATURE_ID', inplace=True)
group_map = {col: col[:6] for col in data.columns}

# 根据分组映射进行分组，收集每组的列名
from collections import defaultdict
groups = defaultdict(list)
for col, prefix in group_map.items():
    groups[prefix].append(col)
for prefix, cols in groups.items():
    # 根据分组获取相应列的 DataFrame
    group_df = data[cols]
    # 创建文件名
    filename = f'{prefix}.txt'

    group_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/expression/{species}/{filename}', sep='\t')





