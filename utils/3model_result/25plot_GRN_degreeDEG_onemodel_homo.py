import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import sys
import os
import omicverse as ov
from sklearn.preprocessing import MinMaxScaler
from matplotlib.font_manager import FontProperties
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
from icecream import install,ic
install()
# 创建一个2x2的子图布局
#fig= plt.figure(figsize=(10, 10))

data_types=['TF','PPI','KEGG']

exps=['CK','CD','CE']
n=0
# colors =['#556B2F','#556B2F','#556B2F',
#         '#800000', '#800000','#800000',
#         '#191970','#191970','#191970',]
colors = [#(77,187,213), (77,187,213),(77,187,213),
          (0,160,135),(230,75,53),
          (0,160,135),(230,75,53),
          (0,160,135),(230,75,53),
          ]
colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
colors= colors_normalized
gos=['all']#,'nitrogen','photosynthesis','stoma''transpiration',,'water_deprivation'

#for go in gos:
    #n=0
# fig, axs = plt.subplots(3, 3, figsize=(15, 12))
# ax=axs.flatten()
for data_type in data_types:
    # if data_type== 'KEGG':
    #     species=['indica','zeamays','japonica']
    # else:
    #     species=['indica_BGI','zeamays','japonica']
    species=['homo']
    for speci in species:
        for go in gos:
            all_deg_datas=[]
            for exp in exps:
            
    # csv_files = [f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_num_data.txt',
    #              f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_num_transpiration.csv', 
    #              f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_num_nitrogen.csv', 
    #             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_num_photosynthesis.csv',
    #             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_num_water_deprivation.csv']
                csv_files = [f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{exp}_c_degree_{go}.csv',
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_c_closeness_all.csv',
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_pagerank_all.csv',
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_between_all.csv',
                            ]
                dataframes = []
                for file in csv_files:
                    df = pd.read_csv(file)
                    # 只保留第二列数据，假设列名为'value'
                    df = df.iloc[:, [1]]  # 根据实际情况可能需要调整列的选择
                    dataframes.append(df)

                # 假设第一个文件的第一列可以作为合并的基准
                base_df = pd.read_csv(csv_files[0]).iloc[:, [0]]
                combined_df = pd.concat([base_df] + dataframes, axis=1)
                scaler = MinMaxScaler(feature_range=(0, 1))
                # 对除了第一列（基准列）以外的所有列应用缩放
                combined_df.iloc[:, 1:] = scaler.fit_transform(combined_df.iloc[:, 1:])

                # 计算平均值
                combined_df[f'avg_{exp}'] =combined_df.iloc[:, 1] #combined_df.iloc[:, 1:].mean(axis=1)
                all_deg_datas.append(combined_df)
            common_keys = pd.merge(all_deg_datas[0]['Unnamed: 0'], all_deg_datas[1]['Unnamed: 0'], how='inner', on='Unnamed: 0')
            common_keys = pd.merge(common_keys, all_deg_datas[2]['Unnamed: 0'], how='inner', on='Unnamed: 0')
            #ic(common_keys)
            all_deg_datas[0] = all_deg_datas[0][all_deg_datas[0]['Unnamed: 0'].isin(common_keys['Unnamed: 0'])].sort_values('Unnamed: 0')
            all_deg_datas[1] = all_deg_datas[1][all_deg_datas[1]['Unnamed: 0'].isin(common_keys['Unnamed: 0'])].sort_values('Unnamed: 0')
            all_deg_datas[2] = all_deg_datas[2][all_deg_datas[2]['Unnamed: 0'].isin(common_keys['Unnamed: 0'])].sort_values('Unnamed: 0')

            #ic(all_deg_datas)
            all_deg_datas[0]['avg_CK'] = (all_deg_datas[0]['avg_CK'] * 10000).astype(np.float16).tolist()
            all_deg_datas[1]['avg_CD'] = (all_deg_datas[1]['avg_CD'] * 10000).astype(np.float16).tolist()
            all_deg_datas[2]['avg_CE'] = (all_deg_datas[2]['avg_CE'] * 10000).astype(np.float16).tolist()
            ic(all_deg_datas[0]['avg_CK'],all_deg_datas[1]['avg_CD'],all_deg_datas[2]['avg_CE'])

            DEG_data=np.stack([all_deg_datas[0]['avg_CK'],all_deg_datas[1]['avg_CD'],all_deg_datas[2]['avg_CE']],axis=1)
            #DEG_data=np.stack([all_deg_datas[0]['avg_CK'],all_deg_datas[1]['avg_CD']],axis=1)

            #DEG_data=DEG_data.fillna(53)
            #ic(DEG_data)
            #DEG_data['avg_CE'] = DEG_data['avg_CE'].astype(int)
            DEG_data=pd.DataFrame(DEG_data,columns=['avg_CK','avg_CD','avg_CE'],index=common_keys['Unnamed: 0'])
            #DEG_data=pd.DataFrame(DEG_data,columns=['avg_CK','avg_CD'],index=common_keys['Unnamed: 0'])
            #DEG_data=DEG_data.fillna(4000)
            ic(DEG_data)
            for col in list(DEG_data.columns):  # 创建一个列名列表的副本，因为我们将改变 DataFrame
                for i in range(1, 4):  # 从1到3，因为我们想要总共3个副本（包括原始的）
                    DEG_data[f'{col}.{i}'] = DEG_data[col]
            #ic(DEG_data)
            DEG_data.iloc[:,1:] = DEG_data.iloc[:,1:].astype('float') + 1e-16
            ck_columns = DEG_data.filter(like='CK', axis=1).reset_index(drop=True)#pd.DataFrame(DEG_data['avg_CK'],columns=['avg_CK'])
            cd_columns = DEG_data.filter(like='CD', axis=1).reset_index(drop=True)#pd.DataFrame(DEG_data['avg_CD'],columns=['avg_CD'])
            ce_columns = DEG_data.filter(like='CE', axis=1).reset_index(drop=True)#pd.DataFrame(DEG_data['avg_CE'],columns=['avg_CE'])
            #ic(ck_columns)
            deg_data_ck_cd =  pd.concat([ck_columns, cd_columns], axis=1).astype(float)
            deg_data_ck_ce =  pd.concat([ck_columns, ce_columns], axis=1).astype(float)
            deg_data_cd_ce = pd.concat([cd_columns, ce_columns], axis=1).astype(float)

            dds_ck_cd=ov.bulk.pyDEG(deg_data_ck_cd)
            dds_ck_cd.drop_duplicates_index()
            print('... drop_duplicates_index success')
            dds_ck_cd.normalize()
            print('... estimateSizeFactors and normalize success')

            dds_ck_ce=ov.bulk.pyDEG(deg_data_ck_ce)
            dds_ck_ce.drop_duplicates_index()
            print('... drop_duplicates_index success')
            dds_ck_ce.normalize()
            print('... estimateSizeFactors and normalize success')

            dds_cd_ce=ov.bulk.pyDEG(deg_data_cd_ce)
            dds_cd_ce.drop_duplicates_index()
            print('... drop_duplicates_index success')
            dds_cd_ce.normalize()
            print('... estimateSizeFactors and normalize success')

            # cd_columns = cd_columns.astype('float16')
            # ce_columns = ce_columns.astype('float16')
            # ck_columns = ck_columns.astype('float16')
            # ic(cd_columns,ck_columns,ce_columns)
            treatment_groups_cd=pd.Index(['avg_CD','avg_CD.1','avg_CD.2','avg_CD.3'])#cd_columns#.columns#.tolist()
            treatment_groups_ce=pd.Index(['avg_CE','avg_CE.1','avg_CE.2','avg_CE.3'])#ce_columns#.columns#.tolist()
            control_groups_ck=pd.Index(['avg_CK','avg_CK.1','avg_CK.2','avg_CK.3'])#ck_columns#.columns#.tolist()
            #ic(treatment_groups_cd,treatment_groups_ce,control_groups_ck)
            ic(dds_ck_cd)
            result_ck_cd=dds_ck_cd.deg_analysis(treatment_groups_cd,control_groups_ck,method='ttest')
            #ic(result_ck_cd)
            result_ck_cd = pd.DataFrame(result_ck_cd)
            combined_df_ck_cd = pd.concat([common_keys['Unnamed: 0'], result_ck_cd], axis=1)

            result_ck_ce=dds_ck_ce.deg_analysis(treatment_groups_ce,control_groups_ck,method='ttest')
            result_ck_ce = pd.DataFrame(result_ck_ce)
            combined_df_ck_ce = pd.concat([common_keys['Unnamed: 0'], result_ck_ce], axis=1)

            result_cd_ce=dds_cd_ce.deg_analysis(treatment_groups_cd,treatment_groups_ce,method='ttest')
            result_cd_ce = pd.DataFrame(result_cd_ce)
            combined_df_cd_ce = pd.concat([common_keys['Unnamed: 0'], result_cd_ce], axis=1)

            def label_data(row):
                if row['sig'] == 'sig' and row['log2FC'] > 0:
                    return 1
                elif row['sig'] == 'sig' and row['log2FC'] < 0:
                    return -1
                else:
                    return 0

            # 在 df 上应用该函数，并将结果赋给 'Label' 列
            # combined_df_ck_cd['Label'] = combined_df_ck_cd.apply(label_data, axis=1)
            # if test==1:
            #     exp_data = pd.read_csv(f'{exp_datapath}')
            # else:
            #     with open(f'{exp_datapath}','rb') as fp:
            #         exp_data=pickle.load(fp)

            #exp_data['Label'] =combined_df_ck_cd.apply(label_data, axis=1)
            directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CD/'
            if not os.path.exists(directory):
                os.makedirs(directory)
            combined_df_ck_cd.to_csv(f'{directory}/{go}_{speci}_CK_CD_DEG_alloutput.csv',index=False)
            ic(combined_df_ck_cd)
            df = combined_df_ck_cd
            #df['-log10pvalue'] = -np.log10(df['pvalue'])
            # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
            #selected_genes = df[(df['-log10pvalue'] > 0.05) & (df['log2FC'].abs() > 0.585)]
            # selected_genes_up = df[(df['pvalue'] < 0.05) & (df['abs(log2FC)'] > 0.585)]  #1.5倍
            # sorted_selected_genes_up = selected_genes_up.sort_values(by='abs(log2FC)', ascending=False)
            # sorted_selected_genes_up.to_csv(f'{directory}/{go}_{speci}_CK_CD_DEG.csv',index=False)

            selected_genes_up = df[(df['qvalue'] < 0.05) & (df['log2FC'] > 0.585)]
            selected_genes_up.to_csv(f'{directory}/{speci}_CK_CD_DEG_up.csv',index=False)
            selected_genes_down = df[(df['qvalue'] < 0.05) & (df['log2FC'] < -0.585)]
            selected_genes_down.to_csv(f'{directory}/{speci}_CK_CD_DEG_down.csv',index=False)

            all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
            sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
            sorted_selected_genes_up.to_csv(f'{directory}/{go}_{speci}_CK_CD_DEG.csv',index=False)

            #打印结果
            print(f"{speci}CK_CD上调基因数量：{len(selected_genes_up)}")
            print(f"{speci}CK_CD下调基因数量：{len(selected_genes_down)}")

            directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CE/'
            if not os.path.exists(directory):
                os.makedirs(directory)
            combined_df_ck_ce['Label'] = combined_df_ck_ce.apply(label_data, axis=1)
            # exp_data = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv')
            # exp_data['Label'] =combined_df_ck_ce.apply(label_data, axis=1)
            combined_df_ck_ce.to_csv(f'{directory}/{speci}_CK_CE_DEG_alloutput.csv',index=False)

            df = combined_df_ck_ce
            df['-log10pvalue'] = -np.log10(df['pvalue'])
            # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
            selected_genes_up = df[(df['qvalue'] < 0.05) & (df['log2FC'] > 0.585)]
            selected_genes_up.to_csv(f'{directory}/{speci}_CK_CE_DEG_up.csv',index=False)
            selected_genes_down = df[(df['qvalue'] < 0.05) & (df['log2FC'] < -0.585)]
            selected_genes_down.to_csv(f'{directory}/{speci}_CK_CE_DEG_down.csv',index=False)
            print(f"{speci}CK_CE上调基因数量：{len(selected_genes_up)}")
            print(f"{speci}CK_CE下调基因数量：{len(selected_genes_down)}")
            all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
            sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
            sorted_selected_genes_up.to_csv(f'{directory}/{go}_{speci}_CK_CE_DEG.csv',index=False)

            directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CD_CE/'
            if not os.path.exists(directory):
                os.makedirs(directory)
            combined_df_cd_ce['Label'] = combined_df_cd_ce.apply(label_data, axis=1)
            # exp_data = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv')
            # exp_data['Label'] =combined_df_ck_ce.apply(label_data, axis=1)
            combined_df_cd_ce.to_csv(f'{directory}/{speci}_CD_CE_DEG_alloutput.csv',index=False)

            df = combined_df_cd_ce
            df['-log10pvalue'] = -np.log10(df['pvalue'])
            # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
            #selected_genes = df[(df['-log10pvalue'] > 0.05) & (df['log2FC'].abs() > 0.585)]
            # 根据abs(log2FC)降序排列
            selected_genes_up = df[(df['qvalue'] < 0.05) & (df['log2FC'] > 0.585)]
            selected_genes_up.to_csv(f'{directory}/{speci}_CD_CE_DEG_up.csv',index=False)
            selected_genes_down = df[(df['qvalue'] < 0.05) & (df['log2FC'] < -0.585)]
            selected_genes_down.to_csv(f'{directory}/{speci}_CD_CE_DEG_down.csv',index=False)
            print(f"{speci}CD_CE上调基因数量：{len(selected_genes_up)}")
            print(f"{speci}CD_CE下调基因数量：{len(selected_genes_down)}")
            all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
            sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
            sorted_selected_genes_up.to_csv(f'{directory}/{go}_{speci}_CD_CE_DEG.csv',index=False)






# title_fontsize = 30
# axis_label_fontsize = 25
# title_font = FontProperties()
# title_font.set_weight('bold')
# def calculate_diff_df(df,col1,col2):
#     df['diff_times'] = (df[[col1, col2]].max(axis=1) + 1e-9) / (df[[col1, col2]].min(axis=1) + 1e-9)
#     df_sorted = df.sort_values(by='diff_times', ascending=False)
#     df_sorted = df_sorted.drop(columns=['diff_times'])
#     return df_sorted
# species= sys.argv[1]
# conditions=['CK','CD','CE']

# interested_GOterm=pd.read_excel('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_enrichr/GO/multi_GO_result_trans.xlsx',sheet_name='interest')
# for GO in interested_GOterm['GO']:

#     all_merged = []
#     GRN_cluster_density=[]
#     for condition in conditions:
#         file_path=f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{condition}_{GO.replace("/", "")}.csv'
#         data = pd.read_csv(file_path)
#         if data.empty:
#             continue
#         G = nx.DiGraph()
#         G=nx.from_pandas_edgelist(data, 'TF', 'Target',edge_attr='pred_weights1',create_using=nx.DiGraph())
#         #PageRank 算法  是特征向量中心性的一个变种
#         pagerank_list = nx.pagerank(G, alpha=0.9)
#         pagerank_df =  pd.DataFrame([pagerank_list]).T
#         #ic(pagerank_df)

#         #计算点度中心性(degree centrality)作为对比  度中心度
#         print("degree centrality")
#         c_degree = nx.degree_centrality(G)
#         c_degree_value = list(c_degree.values())
#         c_degree_key = list(c_degree.keys())
#         c_degree_df=pd.DataFrame([c_degree]).T
#         #ic(c_degree_df)

#         #计算紧密中心度接近度中心度
#         print("Closeness centrality")
#         c_closeness = nx.closeness_centrality(G)
#         c_closeness_value = list(c_closeness.values())
#         c_closeness_key = list(c_closeness.keys())
#         c_closeness_df=pd.DataFrame([c_closeness]).T
#         #ic(c_closeness_df)

#         # if GO=='regulation of response to water deprivation' or 'positive regulation of response to water deprivation':
#         #     katz_centrality = nx.betweenness_centrality(G, k=int(len(G.nodes()) * 0.8))
#         #     katz_centrality_df=pd.DataFrame([katz_centrality]).T
#         # else:
#         katz_centrality = nx.katz_centrality(G)
#         katz_centrality_df=pd.DataFrame([katz_centrality]).T
#         #介数中心度   中间中心度 中介中心度（Betweenness Centrality）
#         # print("Betweenness Centrality")
#         # #k_value = max(int(len(G.nodes()) * 0.1), 1)
#         # #c_betweenness = nx.betweenness_centrality(G, k=k_value, normalized=True)
#         # c_betweenness = nx.betweenness_centrality(G, k=int(len(G.nodes()) * 0.8))
#         # #c_betweenness = nx.betweenness_centrality(G)
#         # c_betweenness_df=pd.DataFrame([c_betweenness]).T
#         # ic(c_betweenness)
#         #计算特征向量中心度(eigenvector centrality)并输出值
#         # print("eigenvector centrality")
#         # c_eigenvector = nx.eigenvector_centrality(G,max_iter=100,tol=0.01)
#         # c_eigenvector_df=pd.DataFrame([c_eigenvector]).T
#         df_merged = pd.concat([pagerank_df,c_degree_df, c_closeness_df,katz_centrality_df], axis=1)
#         df_merged.columns = ['pagerank','degree', 'closeness','katz']
#         #df_merged.columns = 'c_degree'+'c_closeness'
#         #计算网络的聚类系数  可以单独计算比较每个通路的
        
#         #group_centralty = nx.group_betweenness_centrality(G,G.nodes(),normalized=True)
#         #clustering_coefficient = nx.average_clustering(G)
#         # 计算网络的密度
#         density = nx.density(G)
#         #ic(clustering_coefficient,density)
#         #ic(df_merged)
#         #GRN_cluster_density.append(group_centralty)
#         GRN_cluster_density.append(density)
#         all_merged.append(df_merged)

#     if len(all_merged)==0:
#         continue
#     all_data = pd.concat(all_merged, axis=1)
#     if all_data.shape[1] != 12:
#         continue
#     #ic(GRN_cluster_density)
#     ic(all_data.shape)
#     GRN_data = GRN_cluster_density#pd.concat(GRN_cluster_density, axis=1)
#     all_data.columns = ['CK_pagerank','CK_degree', 'CK_closeness','CK_katz',#'CK_c_betweenness','CK_c_eigenvector',
#                         'CD_pagerank','CD_degree', 'CD_closeness','CD_katz',#'CD_c_betweenness','CD_c_eigenvector',
#                         'CE_pagerank','CE_degree', 'CE_closeness','CE_katz',#'CE_c_betweenness','CE_c_eigenvector',
#                         ]
#     # GRN_data.columns = ['CK_clustering','CK_density',
#     #                     'CD_clustering','CD_density',
#     #                     'CE_clustering','CE_density',
#     #                     ]
#     all_data = all_data.dropna()
#     #GRN_data = GRN_data.dropna()
#     #ic(all_data)

#     gene_degree = calculate_diff_df(all_data,'CK_degree','CD_degree')
#     gene_closeness = calculate_diff_df(all_data,'CK_closeness','CD_closeness')
#     gene_pagerank = calculate_diff_df(all_data,'CK_pagerank','CD_pagerank')
#     gene_katz = calculate_diff_df(all_data,'CK_katz','CD_katz')

#     df=all_data
#     gene_df_CE_degree = df[((df['CE_degree'] < df['CK_degree']) & (df['CE_degree'] < df['CD_degree'])) |
#                     ((df['CE_degree'] > df['CK_degree']) & (df['CE_degree'] > df['CD_degree']))]
#     gene_df_CE_closeness = df[((df['CE_closeness'] < df['CK_closeness']) & (df['CE_closeness'] < df['CD_closeness'])) |
#                     ((df['CE_closeness'] > df['CK_closeness']) & (df['CE_closeness'] > df['CD_closeness']))]
#     gene_df_CE_pagerank = df[((df['CE_pagerank'] < df['CK_pagerank']) & (df['CE_pagerank'] < df['CD_pagerank'])) |
#                     ((df['CE_pagerank'] > df['CK_pagerank']) & (df['CE_pagerank'] > df['CD_pagerank']))]
#     gene_df_CE_katz = df[((df['CE_katz'] < df['CK_katz']) & (df['CE_katz'] < df['CD_katz'])) |
#                     ((df['CE_katz'] > df['CK_katz']) & (df['CE_katz'] > df['CD_katz']))]
    
#     fig, axs = plt.subplots(2, 4, figsize=(16, 10))
#     axs = axs.flatten()

#     for i,data in enumerate([gene_katz,gene_degree,gene_closeness,gene_pagerank,GRN_data,gene_df_CE_degree,gene_df_CE_closeness,gene_df_CE_pagerank]):
#         # 柱子宽度
#         #axs[0].set_xlabel('Genes',fontsize=axis_label_fontsize,weight='bold')
#         if i==0:
#             width = 0.35
#             top5_df = data.iloc[:5,:]
#             #ic(top5_df)
#             x = np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width/2, top5_df['CK_katz'], width, label='CK',color="#4C9F70")
#             rects2 = axs[i].bar(x + width/2, top5_df['CD_katz'], width, label='CD',color="#F67280")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene katz',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-katz',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#             #axs[i].set_yticks(x)  # 确保 y 轴刻度正确显示
#             #axs[i].set_yticklabels(top5_df.index)  # 假设索引包含想要显示的标签
#         elif i==1:
#             width = 0.35
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width/2, top5_df['CK_degree'], width, label='CK',color="#4C9F70")
#             rects2 = axs[i].bar(x + width/2, top5_df['CD_degree'], width, label='CD',color="#F67280")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Degrees',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-degree',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#         elif i==2:
#             width = 0.35
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width/2, top5_df['CK_closeness'], width, label='CK',color="#4C9F70")
#             rects2 = axs[i].bar(x + width/2, top5_df['CD_closeness'], width, label='CD',color="#F67280")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Closeness',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-closeness',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#         elif i==3:
#             width = 0.35
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width/2, top5_df['CK_pagerank'], width, label='CK',color="#4C9F70")
#             rects2 = axs[i].bar(x + width/2, top5_df['CD_pagerank'], width, label='CD',color="#F67280")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Pagerank',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-pagerank',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#         elif i==4:
#             width = 0.25
#             top5_df = data#.iloc[:5,:]
#             #ic(top5_df)
#             x = 1#np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width, top5_df[0], width, label='CK', color="#4C9F70")
#             rects2 = axs[i].bar(x, top5_df[1], width, label='CD', color="#F67280")
#             rects3 = axs[i].bar(x + width, top5_df[2], width, label='CE', color="#6A4C93")
#             # 添加一些文本标签和标题
#             axs[i].set_xlabel('GRN density', fontsize=axis_label_fontsize, weight='bold')
#             axs[i].set_title('CK-CD-CE density', fontsize=title_fontsize, weight='bold')
#             #axs[i].set_yticks(x)  # 确保 y 轴刻度正确显示
#             #axs[i].set_yticklabels(top5_df.index)  # 假设索引包含想要显示的标签
#         elif i==5:
#             width = 0.25
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))  # 基因标签位置
#             rects1 = axs[i].bar(x - width, top5_df['CK_degree'], width, label='CK', color="#4C9F70")
#             rects2 = axs[i].bar(x, top5_df['CD_degree'], width, label='CD', color="#F67280")
#             rects3 = axs[i].bar(x + width, top5_df['CE_degree'], width, label='CE', color="#F8B595")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Degrees',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-CE-degree',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#         elif i==6:
#             width = 0.25
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))
#             rects1 = axs[i].bar(x - width, top5_df['CK_closeness'], width, label='CK', color="#4C9F70")
#             rects2 = axs[i].bar(x, top5_df['CD_closeness'], width, label='CD', color="#F67280")
#             rects3 = axs[i].bar(x + width, top5_df['CE_closeness'], width, label='CE', color="#F8B595")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Closeness',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-CE-closeness',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)
#         elif i==7:
#             width = 0.25
#             top5_df = data.iloc[:5,:]
#             x = np.arange(len(top5_df))
#             rects1 = axs[i].bar(x - width, top5_df['CK_pagerank'], width, label='CK', color="#4C9F70")
#             rects2 = axs[i].bar(x, top5_df['CD_pagerank'], width, label='CD', color="#F67280")
#             rects3 = axs[i].bar(x + width, top5_df['CE_pagerank'], width, label='CE', color="#F8B595")
#             # 添加一些文本标签和标题
#             axs[i].set_ylabel('gene Pagerank',fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_title('CK-CD-CE-pagerank',fontsize=title_fontsize,weight='bold')
#             axs[i].set_xticklabels(top5_df.index,rotation=90,fontsize=axis_label_fontsize,weight='bold')
#             axs[i].set_xticks(x)

        

#         #axs[i].legend()
#         axs[i].spines['bottom'].set_linewidth(1.5)  # 设置边框线宽为2.0
#         axs[i].spines['right'].set_linewidth(1.5)
#         axs[i].spines['left'].set_linewidth(1.5)
#         axs[i].spines['top'].set_linewidth(1.5)
#         axs[i].tick_params(labelsize=20,bottom=False, top=False, left=True, right=False)
#         axs[i].legend(fontsize=15, loc='lower right')
#         for label in axs[i].get_xticklabels() + axs[i].get_yticklabels():
#             label.set_weight('bold')
#     ic(GO)
#     plt.title(f'{GO}',fontdict={'fontname': 'Times New Roman', 'fontsize': 16})
#     plt.tight_layout()
#     plt.savefig(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{GO.replace("/", "")}.png',bbox_inches='tight')








