import numpy as np
import pandas as pd
from icecream import install,ic
install()
import sys
import os
from sklearn.metrics import roc_curve, auc
import networkx as nx
# 加载第一个 .npz 文件

#organ = sys.argv[3]
def get_GRN(condition,speci,data_type):
    with np.load(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/train_out.npz') as data:
        labels = data['labels']  # 假设我们使用的关键字是 'arr_0'
        preds1 = data['preds'] 
        # print(labels)
        # print(preds)
        #print('File1 - arr_0:\n', arr1)

    # 加载第二个 .npz 文件
    with np.load(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/val_out.npz') as data:
        labels = data['labels']  # 假设我们使用的关键字是 'arr_0'
        preds2 = data['preds'] 
        #arr2 = data['arr_0']  # 假设我们使用的关键字是 'arr_0'
        #print('File2 - arr_0:\n', arr2)

    # 加载第三个 .npz 文件
    with np.load(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/test_out.npz') as data:
        labels = data['labels']  # 假设我们使用的关键字是 'arr_0'
        preds3 = data['preds']
        #arr3 = data['arr_0']  # 假设我们使用的关键字是 'arr_0'
        #print('File3 - arr_0:\n', arr3)

    # 假设我们想要沿第一个轴（垂直方向）拼接数组
    concatenated_preds = np.vstack((preds1, preds2, preds3))
    #print('Concatenated array:\n', concatenated_array)
    print(concatenated_preds.shape)

    # 载入第一个 CSV 文件，保留标题行
    df1 = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/Train_set.csv')
    df2 = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/Validation_set.csv')
    df3 = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/Test_set.csv')
    ic(df1.shape,df2.shape,df3.shape)
    common_columns = df1.columns.intersection(df2.columns).intersection(df3.columns)

    # 保持列的顺序一致
    df1 = df1[common_columns]
    df2 = df2[common_columns]
    df3 = df3[common_columns]

    # 现在执行纵向拼接
    #concatenated_df = pd.concat([df1, df2, df3], axis=0, ignore_index=True)
    # 将这三个 DataFrame 对象拼接为一个，axis=0 表示沿着行的方向拼接
    concatenated_traindata = pd.concat([df1, df2, df3], axis=0, ignore_index=True)
    ic(concatenated_traindata.shape)

    ic(concatenated_preds)
    if data_type=='PPI':
        preds = pd.DataFrame(np.hstack([concatenated_preds, concatenated_preds]), columns=['TF', 'Target'])
    else:
        preds = pd.DataFrame(concatenated_preds, columns=['TF', 'Target'])
    #print(array_df)
    # 现在，使用 pandas 的 concat 函数将两个 DataFrame 按行拼接起来
    all_data = pd.concat([concatenated_traindata, preds], axis=1, ignore_index=True)
    #print(all_data)

    gene_df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}//ALL_Genefile.csv')  # 替换为你的 CSV 文件路径

    gene_map = pd.Series(gene_df.Gene.values, index=gene_df.index).to_dict()

    #ic(all_data)
    all_data.iloc[:,1] = all_data.iloc[:,1].map(gene_map)
    all_data.iloc[:,2] = all_data.iloc[:,2].map(gene_map)
    # df['第一列'] = df['第一列'].fillna('默认值')
    # df['第二列'] = df['第二列'].fillna('默认值')

    #ic(all_data)

    #----------------------------------生成LOC-OS的TF家族-----------------------------------------------
    if speci == 'zeamays':
        relation={}
        for i in open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt"):
            rap=str(i.split()[0])
            msu=str(i.split()[1])
            #ic(rap)
            #ic(msu)
            if msu!="None":
                if " " in msu:
                    for a in msu.split(" "):
                        relation[a] = rap
                else:
                    relation[msu] = rap
        #ic(relation)
        # 读取目标 CSV 文件中的某一列为 DataFrame
        input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zma_TF_list.txt'  # 替换为真实的 csv 文件路径
        df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
        def invert_dict(original_dict):
            """翻转字典，使得原始的值变成键，原始的键变成值。"""
            inverted_dict = {}
            for key, value in original_dict.items():
                if value in inverted_dict:
                    # 确保字典的值是唯一的，如果不唯一，则添加到列表中
                    if not isinstance(inverted_dict[value], list):
                        inverted_dict[value] = [inverted_dict[value]]
                    inverted_dict[value].append(key)
                else:
                    inverted_dict[value] = key
            return inverted_dict

        def map_values_to_keys(value_list, mapping_dict):
            """使用翻转后的字典映射列表中的值到对应的键。"""
            inverted_dict = invert_dict(mapping_dict)
            return [inverted_dict.get(value, None) for value in value_list]
        data1 = np.array(df.iloc[:, 1].apply(lambda x: f'{x}'))
        #ic(data1)
        mapped_keys = map_values_to_keys(data1,relation)
        #ic(mapped_keys)
        #ic(df.iloc[:,1])
        # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
        #keys_found = [key for key, val in relation.items() if val in df.iloc[:,1]]
        df['Os_geneid'] = mapped_keys#df.iloc[:,1].map(relation).fillna("None")
        #ic(mapped_keys)
        #ic(df['Os_geneid'])
    elif speci=='indica_BGI':

        input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Indica_Osi_TF_list.txt'  # 替换为真实的 csv 文件路径
        df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
        #ic(df)
        # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
        df['Os_geneid'] = df.iloc[:,1]

    else:
        relation={}
        for i in open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data_proceed/MSU2RAP/RAP-MSU_2018-03-29.txt"):
            rap=str(i.split()[0])
            msu=str(i.split()[1])
            if msu!="None":
                if "," in msu:
                    for a in msu.split(","):
                        relation[a[0:-2]] = rap
                else:
                    relation[msu[0:-2]] = rap
        #ic(relation)
        # 读取目标 CSV 文件中的某一列为 DataFrame
        input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/LOC_OS_TF_list.csv'  # 替换为真实的 csv 文件路径
        df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
        #ic(df)
        # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
        df['Os_geneid'] = df.iloc[:,1].map(relation).fillna("None")
        #ic(df['Os_geneid'])
        #ic(len(df['Os_geneid'].unique()))

    # # 定义输出 CSV 文件的路径
    # output_csv_path = '/path/to/your/output.csv'  # 替换为你想要保存 csv 文件的路径

    # # 将结果写入新的 csv 文件
    # df.to_csv(output_csv_path, index=False)

    # 假设 all_data 和 df 已经定义了，并且 all_data 的第一列没有列名，我们先给它命名为 'ID' 方便操作
    all_data.columns =['ID'] +  ['TF'] + ['Target'] +['label'] +['pred_weights1'] +['pred_weights2'] # 如果all_data的第一列已经有列名，跳过这一步
    all_data['TF'] = all_data['TF'].astype(str)
    all_data['Target'] = all_data['Target'].astype(str)
    all_data['ID'] = all_data['ID'].astype(str)
    all_data['label'] = all_data['label'].astype(str)
    all_data['pred_weights1'] = all_data['pred_weights1'].astype(str)
    all_data['pred_weights2'] = all_data['pred_weights2'].astype(str)
    # 用 pd.merge 进行左连接合并，基于 all_data 的 'ID' 和 df 的 'Os_geneid'
    #ic(all_data)
    merged_datatf = pd.merge(all_data, df[['Os_geneid', 'Family']], left_on='TF', right_on='Os_geneid', how='left')
    merged_datatarget = pd.merge(all_data, df[['Os_geneid', 'Family']], left_on='Target', right_on='Os_geneid', how='left')
    #ic(merged_datatf)
    counter = 0  # 自增计数器
    def generate_target_gene():
        #global counter
        #counter += 1{counter}
        return f'target_gene'

    # Apply 函数来填充空值
    merged_datatf['Family'] = merged_datatf['Family'].apply(lambda x: generate_target_gene() if pd.isna(x) else x)
    merged_datatarget['Family'] = merged_datatarget['Family'].apply(lambda x: generate_target_gene() if pd.isna(x) else x)
    #ic(merged_datatf)
    # 将不匹配的值（NaN值）替换为 'target_gene'
    #merged_data['Family'].fillna('target_gene', inplace=True)

    # 现在 'Family' 列包含了匹配的 Family 值 或者 'target_gene'
    # 如果不需要保留 'Os_geneid' 列，可以将它从结果中删除
    merged_datatf.drop(columns=['Os_geneid'], inplace=True)
    merged_datatarget.drop(columns=['Os_geneid'], inplace=True)
    #ic(merged_datatf)
    # 如果想要只保留 'ID' 和新的 'Family' 列，你可以进一步选择这些特定的列：
    resulttf = merged_datatf[['TF', 'Family']]
    resulttarget = merged_datatarget[['Target', 'Family']]
    #ic(resulttf)
    #ic(all_data)

    tf_to_family = resulttf.set_index('TF')['Family'].to_dict()
    all_data['tf_family'] = all_data['TF'].map(tf_to_family)
    target_to_family = resulttarget.set_index('Target')['Family'].to_dict()
    all_data['target_family'] = all_data['Target'].map(target_to_family)


    # all_data = pd.merge(all_data, resulttf[['TF', 'Family']], left_on='TF', right_on='TF', how='left')#resulttf['Family']
    # all_data = pd.merge(all_data, resulttarget[['Target', 'Family']], left_on='Target', right_on='Target', how='left')#resulttarget['Family']
    #['target_family']['tf_family']
    #ic(all_data['family'].unique())
    ic(all_data)
    combined_series = pd.concat([all_data['TF'].astype(str), all_data['Target'].astype(str)])
    unique_values = combined_series.unique()
    #ic(len(unique_values))
    #-------------------------添加差异表达信息---------------------------
    # # 读取 DEG_output.csv 文件
    # DEG_output_ck_cd = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/{species}_CK_CD_DEG_output.csv')
    # DEG_output_ck_ce = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/{species}_CK_CE_DEG_output.csv')
    # # 将 DEG_output 的 'GeneID' 列设为索引，便于查找
    # DEG_output_ck_cd.set_index('FEATURE_ID', inplace=True)
    # DEG_output_ck_ce.set_index('FEATURE_ID', inplace=True)
    # # 使用 map 函数基于 'Target' 列匹配 'Label' 值
    # # 这里假设 'Target' 列在 all_data 中，'Label' 在 DEG_output 中
    # ic(all_data,DEG_output_ck_cd)
    # ic(DEG_output_ck_cd['Label'].index.is_unique)
    # DEG_output_ck_cd_unique_label = DEG_output_ck_cd['Label'][~DEG_output_ck_cd['Label'].index.duplicated(keep='first')]
    # DEG_output_ck_cd_unique_log2fc = DEG_output_ck_cd['log2FC'][~DEG_output_ck_cd['log2FC'].index.duplicated(keep='first')]

    # # DEG_output_ck_cd_unique = DEG_output_ck_cd['Label'][~DEG_output_ck_cd['Label'].index.duplicated(keep='first')]
    # # all_data['Target_DEG_ck_cd'] = all_data['Target'].map(DEG_output_ck_cd_unique)

    # # DEG_output_ck_cd_unique = DEG_output_ck_cd['log2FC'][~DEG_output_ck_cd['log2FC'].index.duplicated(keep='first')]
    # # all_data['Target_DEG_log2FC_ck_cd'] = all_data['Target'].map(DEG_output_ck_cd_unique)

    # all_data['Target_DEG_ck_cd'] = all_data['Target'].map(DEG_output_ck_cd_unique_label)
    # all_data['Target_DEG_log2FC_ck_cd'] = all_data['Target'].map(DEG_output_ck_cd_unique_log2fc)

    # all_data['TF_DEG_ck_cd'] = all_data['TF'].map(DEG_output_ck_cd_unique_label)
    # all_data['TF_DEG_log2FC_ck_cd'] = all_data['TF'].map(DEG_output_ck_cd_unique_log2fc)

    # all_data['Target_DEG_ck_ce'] = all_data['Target'].map(DEG_output_ck_cd_unique_label)
    # all_data['Target_DEG_log2FC_ck_ce'] = all_data['Target'].map(DEG_output_ck_cd_unique_log2fc)

    # all_data['TF_DEG_ck_ce'] = all_data['TF'].map(DEG_output_ck_cd_unique_label)
    # all_data['TF_DEG_log2FC_ck_ce'] = all_data['TF'].map(DEG_output_ck_cd_unique_log2fc)

    #ic(all_data)
    #-------------------------添加GO富集信息---------------------------
    # 读取 OSGO.csv 文件
    if speci =='indica':
        OSGO = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/japonica/OSGO.csv')
    else:
        OSGO = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/{speci}/OSGO.csv')
    # 将 OSGO 数据按 OSGeneid 分组，然后将每个组的 GO_ID 聚合成列表
    GO_ID_map = OSGO.groupby('OSGeneid')['GO_ID'].apply(list).to_dict()
    GO_term_map = OSGO.groupby('OSGeneid')['GO_term'].apply(list).to_dict()

    # 用映射填充 all_data 中的 Target_GO 列
    # 如果 Target 在 GO_ID_map 中找到匹配，就填充对应的 GO_ID 列表
    # 如果没有找到匹配，就填充 NaN
    all_data['Target_GO'] = all_data['Target'].map(GO_ID_map)
    all_data['Target_GO'] = all_data['Target_GO'].apply(lambda x: x if isinstance(x, list) else [])

    all_data['TF_GO'] = all_data['TF'].map(GO_ID_map)
    all_data['TF_GO'] = all_data['TF_GO'].apply(lambda x: x if isinstance(x, list) else [])

    all_data['Target_term'] = all_data['Target'].map(GO_term_map)
    all_data['Target_term'] = all_data['Target_term'].apply(lambda x: x if isinstance(x, list) else [])

    all_data['TF_term'] = all_data['TF'].map(GO_term_map)
    all_data['TF_term'] = all_data['TF_term'].apply(lambda x: x if isinstance(x, list) else [])

    # 结果展示
    #ic(all_data['TF_term'])
    #---------------------------------添加点的描述信息--------------------------------------------
    # def set_character(val): 
    #     if val.startswith('target_gene'): #and any(char.isdigit() for char in val)
    #         return 'ellipse'
    #     else:
    #         return 'triangle'

    # # 应用逻辑到tf_family列
    # all_data['tfCharacter'] = all_data['tf_family'].apply(set_character)
    # all_data['targetCharacter'] = all_data['target_family'].apply(set_character)

    # # 创建TFGenecolor列的逻辑
    # def set_color(val):
    #     if val == 1:
    #         return 'green'
    #     elif val == -1:
    #         return 'red'
    #     elif val == 0:
    #         return 'grey'
    #     else:
    #         return None  # 如果有其他意外的值，你可以选择返回None或者其他颜色字符串

    # # 应用逻辑到TF_DEG列
    # all_data['TFGenecolor_ck_cd'] = all_data['TF_DEG_ck_cd'].apply(set_color)
    # all_data['TargetGenecolor_ck_cd'] = all_data['Target_DEG_ck_cd'].apply(set_color)
    # all_data['TFGenecolor_ck_ce'] = all_data['TF_DEG_ck_ce'].apply(set_color)
    # all_data['TargetGenecolor_ck_ce'] = all_data['Target_DEG_ck_ce'].apply(set_color)

    #-----------------------------------筛选符合要求的数据-------------------------------------------
    ic(all_data['pred_weights2'],ic(all_data['label']))
    all_data['pred_weights2'] = all_data['pred_weights2'].astype(float)
    all_data['label'] = all_data['label'].astype(int)
    fpr, tpr, thresholds = roc_curve(list(all_data['label']), list(all_data['pred_weights2']))
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    ic(optimal_threshold)
    data = all_data[all_data['pred_weights2'] >= optimal_threshold]#pd.read_csv(f'/mnt/h/study/deep_learning/gene/project/plot/network_geneid_{pathway}.csv')
    #ic(data['TF_term'])
    directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/'
    if not os.path.exists(directory):
        os.makedirs(directory)
    #data.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{condition}_all.txt',sep='\t',index=False)

    def contains_substr(term_list, substring):
        return any(substring in term for term in term_list if isinstance(term_list, list))

    # interested_GOterm=pd.read_excel('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_enrichr/GO/multi_GO_result_trans.xlsx',sheet_name='interest')
    # for GO in interested_GOterm['GO']:
    #     ic(GO)
    #     GO_data = data[data['TF_term'].apply(lambda x: GO in x) |
    #                        data['Target_term'].apply(lambda x: GO in x)
    #                        ]
    #     ic(len(GO_data))
    #     # GO_data = data[
    #     # data['TF_term'].apply(contains_substr, substring=GO) |
    #     # data['Target_term'].apply(contains_substr, substring=GO)]
    #     GO_data.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_{GO.replace("/", "")}.csv', index=False)
# def contains_substr(term, substring):
#     return substring in term
# ic(data)
# photosynthesis_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'photosynthesis')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'photosynthetic')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'photosystem')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'carbon')) |   
#     data['GO_term'].apply(lambda term: contains_substr(term, 'carbohydrate'))
#     ]

# nitrogen_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'nitrogen')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'nitrate')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'ammonia')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'ammonium')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'glutamine'))
# ]

# chloroplast_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'chlorophyll')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'chloroplast'))
# ]

# carbon_nitrogen_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'carbon-nitrogen'))
# ]

# water_deprivation_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'water deprivation'))
# ]

# transpiration_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'transpiration'))
# ]

# proline_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'proline'))
# ]

# stoma_rows = data[
#     data['GO_term'].apply(lambda term: contains_substr(term, 'stoma')) |
#     data['GO_term'].apply(lambda term: contains_substr(term, 'stomatal'))
# ]
    photosynthesis_rows = data[
        data['TF_term'].apply(contains_substr, substring='photosynthesis') |
        data['Target_term'].apply(contains_substr, substring='photosynthesis') |
        data['TF_term'].apply(contains_substr, substring='photosynthetic') |
        data['Target_term'].apply(contains_substr, substring='photosynthetic') |        
        data['TF_term'].apply(contains_substr, substring='photosystem') |
        data['Target_term'].apply(contains_substr, substring='photosystem') |        
        data['TF_term'].apply(contains_substr, substring='carbon') |
        data['Target_term'].apply(contains_substr, substring='carbon') |        
        data['TF_term'].apply(contains_substr, substring='carbohydrate') |
        data['Target_term'].apply(contains_substr, substring='carbohydrate')  
        ]

    nitrogen_rows = data[
        data['TF_term'].apply(contains_substr, substring='nitrogen') |
        data['Target_term'].apply(contains_substr, substring='nitrogen') |
        data['TF_term'].apply(contains_substr, substring='nitrate') |
        data['Target_term'].apply(contains_substr, substring='nitrate') |        
        data['TF_term'].apply(contains_substr, substring='ammonia') |
        data['Target_term'].apply(contains_substr, substring='ammonia') |        
        data['TF_term'].apply(contains_substr, substring='ammonium') |
        data['Target_term'].apply(contains_substr, substring='ammonium') |        
        data['TF_term'].apply(contains_substr, substring='glutamine') |
        data['Target_term'].apply(contains_substr, substring='glutamine')
        ]

    carbon_nitrogen_rows = data[
        data['TF_term'].apply(contains_substr, substring='carbon-nitrogen') |
        data['Target_term'].apply(contains_substr, substring='carbon-nitrogen')]

    water_deprivation_rows = data[
        data['TF_term'].apply(contains_substr, substring='water deprivation') |
        data['Target_term'].apply(contains_substr, substring='water deprivation')]

    transpiration_rows = data[
        data['TF_term'].apply(contains_substr, substring='transpiration') |
        data['Target_term'].apply(contains_substr, substring='transpiration')]

    chloroplast_rows = data[
        data['TF_term'].apply(contains_substr, substring='chloroplast') |
        data['Target_term'].apply(contains_substr, substring='chloroplast')]

    stoma_rows = data[
        data['TF_term'].apply(contains_substr, substring='stoma') |
        data['Target_term'].apply(contains_substr, substring='stoma')]
    #------------------------------------导出gene txt文件------------------------------------------------
    datalist=[data]#,photosynthesis_rows,nitrogen_rows,water_deprivation_rows,transpiration_rows,carbon_nitrogen_rows,chloroplast_rows
    # datanamelist=['all']#,'photosynthesis','nitrogen','carbon_nitrogen','water_deprivation','transpiration','chloroplast'
    # photosynthesis_rows.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/organ/mini_GO/{organ}_{species}_{exp}_gene_photosynthesis_rows.txt', columns=['TF', 'Target', 'pred_weights2'], sep='\t', index=False)
    # nitrogen_rows.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/organ/mini_GO/{organ}_{species}_{exp}_gene_nitrogen_rows.txt', columns=['TF', 'Target', 'pred_weights2'], sep='\t', index=False)
    # transpiration_rows.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/organ/mini_GO/{organ}_{species}_{exp}_gene_transpiration_rows.txt', columns=['TF', 'Target', 'pred_weights2'], sep='\t', index=False)
    # water_deprivation_rows.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/organ/mini_GO/{organ}_{species}_{exp}_gene_water_deprivation_rows.txt', columns=['TF', 'Target', 'pred_weights2'], sep='\t', index=False)
    path='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/'
    photosynthesis_rows.to_csv(f'{path}/{data_type}_{speci}_{condition}_gene_photosynthesis_rows.csv',index=False)
    nitrogen_rows.to_csv(f'{path}/{data_type}_{speci}_{condition}_gene_nitrogen_rows.csv',index=False)
    transpiration_rows.to_csv(f'{path}/{data_type}_{speci}_{condition}_gene_transpiration_rows.csv',index=False)
    water_deprivation_rows.to_csv(f'{path}/{data_type}_{speci}_{condition}_gene_water_deprivation_rows.csv',index=False)
    stoma_rows.to_csv(f'{path}/{data_type}_{speci}_{condition}_gene_stoma_rows.csv',index=False)


    for i,df in enumerate(datalist):
        
    #     # df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_{datanamelist[i]}.txt', sep='\t', index=False)
    #     # tf_df = df[['TF', 'tfCharacter','TFGenecolor_ck_cd','tf_family']]
    #     # # 从原始DataFrame中选取Target和targetCharacter两列
    #     # target_df = df[['Target', 'targetCharacter','TargetGenecolor_ck_cd','target_family']]
    #     # # 重命名列以便拼接后区分
    #     # tf_df.columns = ['GeneID', 'Character','Genecolor_ck_cd','family']
    #     # target_df.columns = ['GeneID', 'Character','Genecolor_ck_cd','family']
    #     # # 在列方向上拼接两个DataFrame
    #     # concatenated_df = pd.concat([tf_df, target_df], axis=0)
    #     # # 移除重复行，获取唯一值
    #     # unique_rows_df = concatenated_df.drop_duplicates().reset_index(drop=True)
    #     # # 显示结果
    #     # ic(unique_rows_df.shape)
    #     # unique_rows_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_description_{datanamelist[i]}.txt',sep='\t', index=False)

    #     #----------------------------------得到Gene以及TF家族的数量-----------------------------------------------
        gene_series = pd.concat([df['TF'], df['Target']]).drop_duplicates(keep='first')
        gene_counts = pd.concat([df['TF'], df['Target']]).value_counts()

        # 创建一个包含GeneID和Gene_num的新DataFrame
        gene_df = pd.DataFrame({'GeneID': gene_series, 'Gene_num': gene_series.map(gene_counts)}).reset_index(drop=True)

        # 合并'tf_family'和'target_family'列，并计算唯一值及其出现次数
        family_series = pd.concat([df['tf_family'], df['target_family']]).drop_duplicates(keep='first')
        family_counts = pd.concat([df['tf_family'], df['target_family']]).value_counts()

        # 创建一个包含familyID和family_num的新DataFrame
        family_df = pd.DataFrame({'familyID': family_series, 'family_num': family_series.map(family_counts)}).reset_index(drop=True)
        gene_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{condition}_gene_num_all.csv',index=False)#{datalist[i]}
        family_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{condition}_family_num_all.csv',index=False)#{datalist[i]}

        G = nx.DiGraph()
        G=nx.from_pandas_edgelist(df, 'TF', 'Target',edge_attr='pred_weights2',create_using=nx.DiGraph())

        #计算点度中心性(degree centrality)作为对比  度中心度  一个人的社会关系越多，他/她就越重要
        print("degree centrality")
        c_degree = nx.degree_centrality(G)
        c_degree_df=pd.DataFrame([c_degree]).T
        #ic(c_degree_df)

    #     # #中介中心性： 如果一个成员处于其他成员的多条最短路径上，那么该成员就是核心成员。
    #     # bet_cen = nx.betweenness_centrality(G,k=10)
    #     # bet_cen_df=pd.DataFrame([bet_cen]).T

    #     # #计算紧密中心度接近度中心度  接近中心性：一个人跟所有其他成员的距离越近，他/她就越重要。
    #     # print("Closeness centrality")
    #     # c_closeness = nx.closeness_centrality(G)
    #     # c_closeness_df=pd.DataFrame([c_closeness]).T
    #     # ic(c_closeness_df)

    #     # #PageRank 算法  是特征向量中心性的一个变种
    #     # print("pagerank")
    #     # pagerank_list = nx.pagerank(G, alpha=0.9)
    #     # pagerank_df =  pd.DataFrame([pagerank_list]).T
    #     # ic(pagerank_df)

    #     #----------特征向量中心性---------------与你连接的人社会关系越多，你就越重要。
    #     # print("eigenvector_centrality")
    #     # eig_cen = nx.eigenvector_centrality(G,max_iter=1000)
    #     # #katz_centrality_df=eig_cen
    #     # # print("katz_centrality")
    #     # # katz_centrality = nx.katz_centrality(G,max_iter=500)
    #     # katz_centrality_df=pd.DataFrame([eig_cen]).T
    #     #--------------------------------------------------------------
    #    # 如果一个网络很集中，那么势必是中心点，中心度高；而边缘点中心度低。如果一个网络很稀疏，那么中心点、边缘点的中心度没有多少差异。网络聚类系数、网络密度。
    #     #计算网络的聚类系数  可以单独计算比较每个通路的
    #     # group_centralty = nx.group_betweenness_centrality(G,G.nodes(),normalized=True)

    #     # print("average_clustering")
    #     # clustering_coefficient = nx.average_clustering(G,weight='pred_weights1')
    #     # # 计算网络的密度
    #     # print("density")
    #     # density = nx.density(G)
    #     # #局部聚类系数量化了其邻居节点相互聚集构成团(完全图)的程度   概括了整个网络(全局)中的集群
    #     # density = nx.density(G) #密度
        #clusterint_coefficient = nx.average_clustering(G) #平均聚集系数/局部聚集系数   

    #     # num_nodes = nx.number_of_nodes(G) #节点数
    #     # num_edges = nx.number_of_edges(G) #连接数
    
    #     # transitivity = nx.transitivity(G) #传递性/全局聚集系数
    #     # reciprocity = nx.reciprocity(G) #互惠性

        c_degree_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{condition}_c_degree_all.csv')
    #     # c_closeness_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{species}/{exp}_c_closeness_{datanamelist[i]}.csv')
    #     # pagerank_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{species}/{exp}_pagerank_{datanamelist[i]}.csv')
    #     # bet_cen_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{species}/{exp}_between_{datanamelist[i]}.csv')
    #     #katz_centrality_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{species}/{exp}_katz_centrality_{datanamelist[i]}.csv')
    #     #grn_feature=np.array([num_nodes,num_edges,transitivity,reciprocity])#np.array([clustering_coefficient,density])

    #     #pd.DataFrame(grn_feature).to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_grn_feature_add_{datanamelist[i]}.csv')

    #     #ic(clustering_coefficient,density)

    #     tf_target_weights = pd.DataFrame({
    #         'GeneID': pd.concat([df['TF'], df['Target']]),
    #         'pred_weights1': pd.concat([df.loc[df['TF'].index, 'pred_weights1'], df.loc[df['Target'].index, 'pred_weights1']])
    #     })
    #     gene_weights_sum = tf_target_weights.groupby('GeneID')['pred_weights1'].sum().reset_index()
    #     # 你可以选择将结果输出到CSV文件
    #     gene_weights_sum.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/organ/{species}/{organ}_{exp}_gene_weights_{datanamelist[i]}.csv', index=False)

        # # 同理，如果你需要对tf_family和target_family进行相似操作
        # family_target_weights = pd.DataFrame({
        #     'familyID': pd.concat([df['tf_family'], df['target_family']]),
        #     'pred_weights1': pd.concat([df.loc[df['tf_family'].index, 'pred_weights1'], df.loc[df['target_family'].index, 'pred_weights1']])
        # })
        # # 计算每个familyID的pred_weights1总和
        # family_weights_sum = family_target_weights.groupby('familyID')['pred_weights1'].sum().reset_index()

        # # 输出结果到CSV文件
        # family_weights_sum.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_{datanamelist[i]}.csv', index=False)

#============================================通路权重变化========================================================
# df['Combined_GO'] = data.apply(lambda row: list(dict.fromkeys(row['Target_GO'] + row['TF_GO'])), axis=1)
# #data.loc[df.index, 'Combined_GO'] = df.apply(lambda row: list(dict.fromkeys(row['Target_GO'] + row['TF_GO'])), axis=1)
# expanded_rows = [(go_id, row['pred_weights2']) for index, row in df.iterrows() for go_id in row['Combined_GO']]
# #ic(expanded_rows)
# # 转换为DataFrame
# expanded_df = pd.DataFrame(expanded_rows, columns=['GOID', 'Weight'])

# # 对GO ID进行分组，并计算每个GO ID的Weight总和
# goid_weights = expanded_df.groupby('GOID')['Weight'].sum().reset_index()
# goid_weights.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{species}/{exp}_goid_weights.csv', index=False)

conditions=['CK','CD','CE']#
data_types=['TF','PPI','KEGG'] #

for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']
    for speci in species:
        for i,condition in enumerate(conditions):
            get_GRN(condition=condition,speci=speci,data_type=data_type)

