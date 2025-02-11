
import pandas as pd
from icecream import install,ic
install()
import ast
# 指定文件夹路径
dir = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/'
target_go_terms = set(['GO:0009399', 'GO:0006808','GO:0071705',
                   'GO:0006807','GO:0019740','GO:0009768',
                   'GO:0010109','GO:0015979','GO:0009767',
                   'GO:0009654','GO:0009522','GO:0009539',
                   'GO:0009523','GO:0009772','GO:2000123',
                   'GO:2000037','GO:0010119','GO:0010375',
                    'GO:1902456','GO:0010118','GO:0010103',
                    'GO:0090333','GO:2000122','GO:0090332',
                    'GO:0010440','GO:1901527','GO:1990069',
                    'GO:0010148','GO:0009414'
                   ])
def get_matching_gos(target_go, tf_go):
    # 将字符串表示的列表安全地解析为列表
    target_go_list = ast.literal_eval(target_go)
    tf_go_list = ast.literal_eval(tf_go)
    # 合并两个列表
    combined_list = target_go_list + tf_go_list
    # 找出所有匹配的 GO 条目，并返回它们
    matching_gos = [go for go in combined_list if go in target_go_terms]
    return ','.join(matching_gos) if matching_gos else None

# 遍历不同的物种和处理类型
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type == 'TF':
        species=['wheat','sorghum','homo_C3','homo_C4']
    else:
        species = ['wheat','sorghum','homo_C3','homo_C4','homo']
    for speci in species:
        for treat in ['CK', 'CD', 'CE']:
            input_file_path = f'{dir}/{data_type}_{speci}/{treat}_all.txt'
            output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_{treat}_new_regulation.csv'
            go_output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_{treat}_go_filtered.csv'
            # 读取数据，假设文件是以制表符分隔的
            df = pd.read_csv(input_file_path, sep='\t')

            # 计算两列之间的差的绝对值并过滤
            df['diff'] = abs(df['label'] - df['pred_weights2'])
            filtered_df = df[df['diff'] > 0.8]

            # 对过滤后的数据按差值大小进行降序排序
            sorted_filtered_df = filtered_df.sort_values(by='diff', ascending=False)
            #ic(sorted_filtered_df['Target_GO'])

            sorted_filtered_df['Matching_GO'] = sorted_filtered_df.apply(
                        lambda row: get_matching_gos(row['Target_GO'], row['TF_GO']), axis=1
                    )
            go_filtered_df = sorted_filtered_df.dropna(subset=['Matching_GO'])

            go_filtered_df.drop(['ID', 'pred_weights1','Target_GO','TF_GO','Target_term','TF_term','diff'], axis=1, inplace=True)
            # 保存过滤并排序后的数据到新文件
            sorted_filtered_df.to_csv(output_file_path, index=False)
            # 保存包含特定 GO 条目的数据到另一个新文件
            go_filtered_df.to_csv(go_output_file_path, index=False)


    for speci in species:
        CK_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_CK_go_filtered.csv')
        CD_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_CD_go_filtered.csv')
        CE_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_CE_go_filtered.csv')
        
        combined_df = pd.concat([CD_data, CK_data])
        CD_data['source'] = 'CD'
        CK_data['source'] = 'CK'

        # 找出所有行，并标记重复行
        combined_df = pd.concat([CD_data, CK_data])
        combined_df['is_duplicate'] = combined_df.duplicated(subset=['TF', 'Target'], keep=False)

        # 筛选出 CD 中独有的行
        unique_to_cd = combined_df[(combined_df['source'] == 'CD') & (combined_df['is_duplicate'] == False)]

        # 移除辅助列
        CD_CK = unique_to_cd.drop(['source', 'is_duplicate'], axis=1)

        ic(CD_CK)
        combined_df = pd.concat([CE_data, CK_data])
        CE_data['source'] = 'CE'
        CK_data['source'] = 'CK'

        # 找出所有行，并标记重复行
        combined_df = pd.concat([CE_data, CK_data])
        combined_df['is_duplicate'] = combined_df.duplicated(subset=['TF', 'Target'], keep=False)

        # 筛选出 CD 中独有的行
        unique_to_cd = combined_df[(combined_df['source'] == 'CE') & (combined_df['is_duplicate'] == False)]

        # 移除辅助列
        CE_CK = unique_to_cd.drop(['source', 'is_duplicate'], axis=1)

        CD_CK.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_CD_CK.csv',index=False)
        CE_CK.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/new_regulation/{data_type}_{speci}_CE_CK.csv',index=False)
















