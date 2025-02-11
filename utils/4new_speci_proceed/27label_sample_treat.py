import pandas as pd
import numpy as np
import pickle
from icecream import install,ic
install()

#excel_df = pd.read_excel("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/data_description.xlsx")
datapath='/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI'
species=['wheat']#,'sorghum','homo','homo_C3','homo_C4' 'japonica','zeamays','indica_BGI''homo','indica'
for speci in species:
    if speci =='zeamays':
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
    elif speci=='homo' or speci=='homo_C3' or speci=='homo_C4':
        excel_df1 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
        excel_df2 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")
        excel_df3 = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_label_data.csv',sep='\t')
        excel_df4 = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_label_data.csv',sep='\t')
        columns_needed = ['SRRid', 'class']
        df1_selected = excel_df1[columns_needed]
        df2_selected = excel_df2[columns_needed]
        df3_selected = excel_df3[columns_needed]
        df4_selected = excel_df4[columns_needed]
        excel_df = pd.concat([df1_selected, df2_selected,df3_selected,df4_selected,], axis=0)
    elif speci=='wheat':
        excel_df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_label_data.csv',sep='\t')
    elif speci=='sorghum':
        excel_df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_label_data.csv',sep='\t')
    else:
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")

    # 读取 CSV 文件
    with open(f'{datapath}/{speci}/ExpressionData_unique_networkfilter.pkl','rb') as fp:
        csv_df = pickle.load(fp)
    #csv_df = pd.read_csv("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv")
    # if speci=='japonica':
    #     columns_to_drop = [col for col in csv_df.columns if col.startswith('SRR247')]
    #     # 删除这些列
    #     csv_df.drop(columns=columns_to_drop, inplace=True)
            
    ic(csv_df)
    # 假设Excel文件中的列名是'SRRid'和'class'
    # 获取SRRid和class列
    srrid_class_df = excel_df[['SRRid', 'class']]
    # 创建SRRid到class的映射字典
    srrid_to_class_map = dict(zip(srrid_class_df['SRRid'], srrid_class_df['class']))
    # 获取CSV文件的列名
    csv_columns = csv_df.columns
    ic(csv_columns)
    # 创建一个空的标签行
    label_row = pd.Series(index=csv_columns, dtype='object')
    # 遍历CSV的列名查找匹配，并分配标签
    for col in csv_columns:
        for srrid in srrid_to_class_map:
            if str(srrid) in col:
                label_row[col] = srrid_to_class_map[srrid]
                break

    # 重新组织DataFrame为新的顺序和数据
    new_csv_df = pd.concat([label_row.to_frame().T, csv_df])#, ignore_index=True
    ic(new_csv_df)
    # # 重置索引
    # new_csv_df.reset_index(drop=True, inplace=True)
    new_csv_df.to_pickle(f'{datapath}/{speci}/label_filter_expression.pkl')
    #new_csv_df.to_csv(f'{datapath}/label_sample_expression.csv', index=False)









