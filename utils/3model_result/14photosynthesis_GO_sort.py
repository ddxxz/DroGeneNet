
import pandas as pd
from collections import Counter
import ast

# data_types=['TF','PPI','KEGG']
# for data_type in data_types:
#     if data_type== 'KEGG':
#         species=['indica','zeamays','japonica']
#     else:
#         species=['indica_BGI','zeamays','japonica']
#     for speci in species:
#         for treat in ['CK','CD','CE']:
#             df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_photosynthesis_rows.csv')
#             df['Target_term'] = df['Target_term'].apply(ast.literal_eval)
#             df['TF_term'] = df['TF_term'].apply(ast.literal_eval)

#             # 将Target_GO和TF_GO按列拼接
#             df['Combined_term'] = df['Target_term'] + df['TF_term']

#             # 展平所有的列表条目
#             all_go_terms = [item for sublist in df['Combined_term'] for item in sublist]

#             # 对条目进行计数
#             counter = Counter(all_go_terms)

#             # 将计数结果按条目数量排序
#             sorted_items = counter.most_common()

#             # 转换为DataFrame
#             sorted_df = pd.DataFrame(sorted_items, columns=['GO_term', 'Count'])
#             sorted_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_photosynthesis_num/{data_type}_{speci}_{treat}_go_sort.csv',index=False)
#             # 显示结果
#             print(sorted_df)


#===================================绘图===========================================
import pandas as pd
import ast
from collections import Counter
import matplotlib.pyplot as plt

data_types = ['TF', 'PPI', 'KEGG']
species_dict = {
    'TF': ['indica_BGI',  'japonica','zeamays',],
    'PPI': ['indica_BGI',  'japonica','zeamays',],
    'KEGG': ['indica',  'japonica','zeamays',]
}

for data_type in data_types:
    species = species_dict[data_type]
    
    fig, ax = plt.subplots(figsize=(12, 8))  # 创建图表
    
    for speci in species:
        treat = 'CD'
        df = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_photosynthesis_rows.csv')
        #deg_up_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_up.csv')
        #deg_down_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_down.csv')
        #deg_gene=pd.concat([deg_up_data['FEATURE_ID'],deg_down_data['FEATURE_ID']],axis=0)
        #df = df[(df['TF'].isin(deg_gene) | df['Target'].isin(deg_gene))]
        
        #df['tf_family'] = df['Target_term'].apply(ast.literal_eval)
        #df['tf_family'] = df['TF_term'].apply(ast.literal_eval)

        # 将Target_GO和TF_GO按列拼接
        all_go_terms = pd.concat([df['tf_family'], df['target_family']],axis=0)


        # 展平所有的列表条目
       # all_go_terms = df['Combined_term']#[item for sublist in df['Combined_term'] for item in sublist]

        # 对条目进行计数
        counter = Counter(all_go_terms)

        # 将计数结果按条目数量排序
        sorted_items = counter.most_common()

        # 转换为DataFrame
        sorted_df = pd.DataFrame(sorted_items, columns=['TF', 'Count'])
        sorted_df = sorted_df.iloc[1:,]
        sorted_df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_photosynthesis_num/{data_type}_{speci}_{treat}_go_sort.csv', index=False)

        # 绘制柱状图
        top_terms = sorted_df.head(10)  # 只绘制前10个GO_term
        ax.bar(top_terms['TF'], top_terms['Count'], label=speci)

    ax.set_title(f'{data_type} - CD')
    ax.set_xlabel('TF family')
    ax.set_ylabel('Count')
    ax.legend()
    plt.xticks(rotation=90)  # X轴标签旋转90度
    plt.tight_layout()
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_photosynthesis_num/{data_type}_CD_go_barplot.png')
    plt.clf()

















