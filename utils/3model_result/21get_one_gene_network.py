
import pandas as pd

treats=['CK','CD','CE']
geneid ='BGIOSGA031434'       #'Os01g0357100'---LOC_Os01g25484    'Os01g0704100'---LOC_Os01g50820   'Os01g0357100'---BGIOSGA003498   'Os01g0704100'---BGIOSGA004285
# 指定文件路径
file_dir = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/TF_indica_BGI/'  # 路径根据你的文件存放位置进行调整
output_dir = '/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/graphormer2Link_alltype/onegene_network/'
for treat in treats:
    # 读取数据
    df = pd.read_csv(f'{file_dir}{treat}_all.txt', sep='\t')  # 假设文件是以制表符分隔的

    # 筛选包含基因名 'Os01g0357100' 的行
    filtered_df = df[(df['TF'] == geneid) | (df['Target'] == geneid)]

    # 保留 'TF', 'Target', 'pred_weights2' 列
    result_df = filtered_df[['TF', 'Target', 'label','pred_weights2']]

    # 将结果保存到新的 CSV 文件
    result_df.to_csv(f'{output_dir}{treat}_{geneid}_net_indica_NRT11B.csv', index=False)
