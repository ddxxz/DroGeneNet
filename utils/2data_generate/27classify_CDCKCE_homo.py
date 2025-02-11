import pandas as pd
import pickle
from icecream import install,ic
install()

for speci in ['homo']: #'homo','indica',
    for data_type in ['TF','PPI','KEGG']:
        datapath = f"/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}"
        # 加载CSV文件
        with open(f'{datapath}/label_filter_expression.pkl','rb') as fp:
            df = pickle.load(fp)

        #df = pd.read_csv('/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/label_sample_expression.csv', header=None)

        # 获取第二行的数据，用于判断列的分类
        second_row = df.iloc[0]
        #ic(second_row)
        # 初始化存储每个分类列索引的字典，包括第一列
        columns_dict = {'CD': ['FEATURE_ID'], 'CK': ['FEATURE_ID'], 'CE': ['FEATURE_ID']}  # 假设第一列的索引是0

        # 遍历第二行，根据值分类列索引，从第二列开始遍历
        for col_index, value in second_row.items():
            if col_index == 0:
                continue  # 跳过第一列
            if value in columns_dict:
                columns_dict[value].append(col_index)
        #ic(columns_dict)
        # 对每个分类，提取对应的列并保存到新的CSV文件中
        for key, col_indexes in columns_dict.items():
            #ic(col_indexes)
            #col_indexes = col_indexes[1:]
            if len(col_indexes) > 1:  # 如果除了第一列之外还有其他列属于这个分类
                # 提取列
                #ic(df)
                filtered_df = df.loc[:, col_indexes]
                
                # 保存到新的CSV文件
                filtered_df = filtered_df.iloc[1:,:]
                filtered_df = filtered_df.dropna(axis=1)
                ic(filtered_df.shape)
                filtered_df.to_csv(f'{datapath}/{key}_expression.csv', index=False)

#print("文件已根据第二行的CD、CK、CE分类保存，并保留了原始的第一列。")
