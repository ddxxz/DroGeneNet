
import pandas as pd
import numpy as np
import pickle
from icecream import install,ic
install()

#excel_df = pd.read_excel("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/data_description.xlsx")
datapath='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out'
species=['homo','homo_C3','homo_C4']##,,'sorghum''homo','japonica','zeamays','indica','indica_BGI'
for speci in species:
    if speci =='zeamays':
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
    elif speci=='homo' or speci=='homo_C3' or speci=='homo_C4':
        excel_df1 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
        excel_df2 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")
        excel_df3 = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Wheat_label_data.csv",sep='\t')
        excel_df4 = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Sorghum_label_data.csv",sep='\t')
        columns_needed = ['SRRid', 'class']
        df1_selected = excel_df1[columns_needed]
        df2_selected = excel_df2[columns_needed]
        df3_selected = excel_df3[columns_needed]
        df4_selected = excel_df4[columns_needed]
        excel_df = pd.concat([df1_selected, df2_selected,df3_selected,df4_selected], axis=0)
    elif speci=='wheat':
        excel_df = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Wheat_label_data.csv",sep='\t')
    elif speci=='sorghum':
        excel_df = pd.read_csv("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/new_speci/Sorghum_label_data.csv",sep='\t')
    else:
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")

    # 读取 CSV 文件
    #=======================================二选一===============================================
    with open(f'{datapath}/expression/{speci}_expression_batchcorrection.pkl','rb') as fp:
        csv_df = pickle.load(fp)

    csv_df = csv_df.dropna(axis=1, how='any')
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

    nan_columns = new_csv_df.columns[new_csv_df.iloc[0].isna()]

    # 删除含有NaN值的列
    if len(nan_columns) > 0:
        new_csv_df.drop(columns=nan_columns, inplace=True)
        #ic("Columns with NaN values in the 0th row have been removed.")
        #ic(new_csv_df)
    else:
        ic("No columns to remove. No NaN values found in the 0th row.")

    ic(new_csv_df)
    # new_csv_df.reset_index(drop=True, inplace=True)
    new_csv_df.to_pickle(f'{datapath}/expression/{speci}_label_sample_expression.pkl')
    #new_csv_df.to_csv(f'{datapath}/label_sample_expression.csv', index=False)

    with open(f'{datapath}/count/{speci}_count_batchcorrection.pkl','rb') as fp:
        csv_df = pickle.load(fp)

    csv_df = csv_df.dropna(axis=1, how='any')
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

    nan_columns = new_csv_df.columns[new_csv_df.iloc[0].isna()]

    # 删除含有NaN值的列
    # if len(nan_columns) > 0:
    #     new_csv_df.drop(columns=nan_columns, inplace=True)
    #     #ic("Columns with NaN values in the 0th row have been removed.")
    #     #ic(new_csv_df)
    # else:
    #     ic("No columns to remove. No NaN values found in the 0th row.")

    #ic(new_csv_df['0'])
    # new_csv_df.reset_index(drop=True, inplace=True)
    new_csv_df.to_pickle(f'{datapath}/count/{speci}_label_sample_count.pkl')









