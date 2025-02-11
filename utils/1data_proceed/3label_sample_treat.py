#---------------------------------------FPKM表达矩阵---------------------------------------
# import pandas as pd
# import numpy as np
# import pickle
# from icecream import install,ic
# install()
# # 读取 Excel 文件
# #excel_df = pd.read_excel("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/data_description.xlsx")
# datapath='/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies_batchcorrection'
# species=['indica_BGI','indica','japonica','zeamays']
# for speci in species:
#     if speci =='zeamays':
#         excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
#     else:
#         excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")

# # 读取 CSV 文件
# with open(f'{datapath}/{speci}/{speci}_expression_batchcorrection.pkl','rb') as fp:
#     csv_df = pickle.load(fp)
# #csv_df = pd.read_csv("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv")
# ic(csv_df)
# # 假设Excel文件中的列名是'SRRid'和'class'
# # 获取SRRid和class列
# srrid_class_df = excel_df[['SRRid', 'class']]
# # 创建SRRid到class的映射字典
# srrid_to_class_map = dict(zip(srrid_class_df['SRRid'], srrid_class_df['class']))
# # 获取CSV文件的列名
# csv_columns = csv_df.columns
# ic(csv_columns)
# # 创建一个空的标签行
# label_row = pd.Series(index=csv_columns, dtype='object')

# # 遍历CSV的列名查找匹配，并分配标签
# for col in csv_columns:
#     for srrid in srrid_to_class_map:
#         if str(srrid) in col:
#             label_row[col] = srrid_to_class_map[srrid]
#             break

# # 重新组织DataFrame为新的顺序和数据
# new_csv_df = pd.concat([label_row.to_frame().T, csv_df])
# ic(new_csv_df.columns)
# # 按标签行的值对列进行排序，这样相同的标签会聚在一起
# # sorted_columns = new_csv_df.iloc[0].sort_values().index
# # new_csv_df = new_csv_df[sorted_columns]
# # ic(new_csv_df.columns)
# # 获取CSV文件的列标题，假设第一列是'GeneName'
# #ic(new_csv_df)
# # csv_columns = csv_df.columns
# # gene_name_column = csv_columns[0] 
# # csv_data_columns = csv_columns[1:]  # 除去'GeneName'列

# #new_csv_df.insert(0, 'Sample', 'Class')
# ic(len(csv_df.columns))
# #ic(new_csv_df)
# if len(csv_df) <= len(new_csv_df) - 1:
#     # 获取 csv_df 的第一列，忽略其索引
#     #ic(csv_df.iloc[:, 0].tolist()[:5])
#     first_column_values = ["FEATURE_ID"] + csv_df.iloc[:, 0].tolist()#[:-2]
#     #ic(first_column_values)
#     #ic(first_column_values[:5])
#     # 将 first_column_values 赋予 another_df 的第一列，从第二行开始
#     new_csv_df.iloc[0:len(first_column_values), 0] = first_column_values
# else:
#     print("csv_df 的第一列长度大于 another_df 第一列长度减去1，无法替换。")


# #new_csv_df.loc[1:, 0] = csv_df[gene_name_column].values
# #new_csv_df = new_csv_df.iloc[:,:-1]
# ic(new_csv_df.columns)
# ic(new_csv_df.iloc[0])

# nan_counts = new_csv_df.iloc[0].isna().sum()
# #print("NaN values count per column:")
# #ic(nan_counts)
# #ic(new_csv_df.iloc[0])

# nan_columns = new_csv_df.columns[new_csv_df.iloc[0].isna()]

# # 删除含有NaN值的列
# if len(nan_columns) > 0:
#     new_csv_df.drop(columns=nan_columns, inplace=True)
#     #ic("Columns with NaN values in the 0th row have been removed.")
#     #ic(new_csv_df)
# else:
#     ic("No columns to remove. No NaN values found in the 0th row.")
# # second_row = new_csv_df.iloc[0]
# # # 查找第二行中的NaN值并替换为"CK"
# # new_csv_df.iloc[1] = second_row.fillna('CK')

# ic(new_csv_df)
# # df = new_csv_df.iloc[1:]
# # new_header = new_csv_df.iloc[0]

# # # 2. 再将这些标题设为DataFrame的列名
# # new_csv_df.columns = new_header

# # # 3. 然后，删除之前作为标题的那一行
# # new_csv_df = new_csv_df.iloc[1:]

# # # 重置索引
# # new_csv_df.reset_index(drop=True, inplace=True)
# new_csv_df.to_pickle(f'{datapath}/label_sample_expression.pkl')
# #new_csv_df.to_csv(f'{datapath}/label_sample_expression.csv', index=False)


#-------------------------------------genecount表达矩阵-----------------------------------------------
import pandas as pd
import numpy as np
import pickle
from icecream import install,ic
install()

#excel_df = pd.read_excel("/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/data_description.xlsx")
datapath='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/'
species=['sorghum']#'japonica','indica_BGI','wheat','zeamays','homo','indica',
for speci in species:
    if speci =='zeamays':
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
    elif speci=='homo':
        excel_df1 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zea_mays_label_data.csv")
        excel_df2 = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")
        columns_needed = ['SRRid', 'class']
        df1_selected = excel_df1[columns_needed]
        df2_selected = excel_df2[columns_needed]
        excel_df = pd.concat([df1_selected, df2_selected], axis=0)
    else:
        excel_df = pd.read_csv("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Oryza_sativa_label_data.csv")

    # 读取 CSV 文件
    with open(f'{datapath}/expression_nofilter/{speci}_expression_batchcorrection.pkl','rb') as fp:
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
    # 按标签行的值对列进行排序，这样相同的标签会聚在一起
    # sorted_columns = new_csv_df.iloc[0].sort_values().index
    # new_csv_df = new_csv_df[sorted_columns]
    # ic(new_csv_df.columns)
    # 获取CSV文件的列标题，假设第一列是'GeneName'
    #ic(new_csv_df)
    # csv_columns = csv_df.columns
    # gene_name_column = csv_columns[0] 
    # csv_data_columns = csv_columns[1:]  # 除去'GeneName'列

    #new_csv_df.insert(0, 'Sample', 'Class')
    #ic(len(csv_df.columns))
    #ic(new_csv_df)
    # if len(csv_df) <= len(new_csv_df) - 1:
    #     # 获取 csv_df 的第一列，忽略其索引
    #     #ic(csv_df.iloc[:, 0].tolist()[:5])
    #     first_column_values =  csv_df.iloc[:, 0].tolist()#[:-2]["FEATURE_ID"] +
    #     #ic(first_column_values)
    #     #ic(first_column_values[:5])
    #     # 将 first_column_values 赋予 another_df 的第一列，从第二行开始
    #     new_csv_df.iloc[0:len(first_column_values), 0] = first_column_values
    # else:
    #     print("csv_df 的第一列长度大于 another_df 第一列长度减去1，无法替换。")


    #new_csv_df.loc[1:, 0] = csv_df[gene_name_column].values
    #new_csv_df = new_csv_df.iloc[:,:-1]
    #ic(new_csv_df.columns)
    #ic(new_csv_df.iloc[0])

    # nan_counts = new_csv_df.iloc[0].isna().sum()
    # #print("NaN values count per column:")
    # #ic(nan_counts)
    # #ic(new_csv_df.iloc[0])

    # nan_columns = new_csv_df.columns[new_csv_df.iloc[0].isna()]

    # # 删除含有NaN值的列
    # if len(nan_columns) > 0:
    #     new_csv_df.drop(columns=nan_columns, inplace=True)
    #     #ic("Columns with NaN values in the 0th row have been removed.")
    #     #ic(new_csv_df)
    # else:
    #     ic("No columns to remove. No NaN values found in the 0th row.")
    # second_row = new_csv_df.iloc[0]
    # # 查找第二行中的NaN值并替换为"CK"
    # new_csv_df.iloc[1] = second_row.fillna('CK')

    ic(new_csv_df)
    # df = new_csv_df.iloc[1:]
    # new_header = new_csv_df.iloc[0]

    # # 2. 再将这些标题设为DataFrame的列名
    # new_csv_df.columns = new_header

    # # 3. 然后，删除之前作为标题的那一行
    # new_csv_df = new_csv_df.iloc[1:]

    # # 重置索引
    # new_csv_df.reset_index(drop=True, inplace=True)
    new_csv_df.to_pickle(f'{datapath}/expression_nofilter/{speci}_label_sample_expression.pkl')
    #new_csv_df.to_csv(f'{datapath}/label_sample_expression.csv', index=False)









