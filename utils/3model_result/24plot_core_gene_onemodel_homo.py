#--------------------------------------------family_num----------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
from icecream import install,ic
install()
import sys
# 创建一个2x2的子图布局


#exp= sys.argv[1]    # 'CE'
# species= ['japonica','indica', 'zeamays']#'homo'sys.argv[1]
# exps=['CK','CD','CE']#


def plot_core_gene(exps,data_types):
    
    for data_type in data_types:
        # if data_type == 'KEGG':
        #     species=['indica','zeamays','japonica']
        # else:
        #     species=['indica_BGI','zeamays','japonica']
        species=['homo']
        if 'homo' in species:
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            colors = [(77,187,213),
                (0,160,135),(230,75,53)
                ]
        elif len(exps)==3:
            fig, axs = plt.subplots(3, 3, figsize=(15, 12))
            colors = [(77,187,213), (77,187,213), (77,187,213),
                (0,160,135),(0,160,135),(0,160,135),
                (230,75,53),(230,75,53),(230,75,53),
                ]

        else:
            fig, axs = plt.subplots(3, 2, figsize=(12, 12))
            # colors = ["#FFA500", "#FFA500","#FFA500", 
            #            "#800080", "#800080","#800080",
            #             "#008000","#008000","#008000",
            #           ]
            colors = [(77,187,213), (77,187,213),
                    (0,160,135),(0,160,135),
                    (230,75,53),(230,75,53),
                    ]
        colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
        colors= colors_normalized
        
        for speci in species:
            n=0
            for exp in exps:
                csv_files = [f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{exp}_family_num_all.csv', 
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_num_transpiration.csv', 
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_num_nitrogen.csv', 
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_num_photosynthesis.csv', 
                            # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_num_water_deprivation.csv'
                            ]

                #title = [f'{exp}-all',f'{exp}-transpiration',f'{exp}-nitrogen',f'{exp}-photosynthesis',f'{exp}-water_deprivation']
                title = f'{data_type}_{exp}-{speci}'

                # 遍历CSV文件和子图轴
                # 定义标题和轴标签的字体大小
                title_fontsize = 25
                axis_label_fontsize = 18
                axs=axs.flatten()

                for i, csv_file in enumerate(csv_files):
                    # 加载CSV文件
                    import pandas as pd
                    ic(csv_file)
                    if csv_file.endswith('.txt'): 
                        df = pd.read_csv(csv_file,sep='\t')
                    else:
                        df = pd.read_csv(csv_file)
                    
                    df = df.sort_values(by='family_num', ascending=False)
                    ic(df)
                    # 删除第一行           
                    df = df.sort_values(by='family_num', ascending=False).iloc[1:7]  # Sort and filter rows
                    # df = df.iloc[1:]
                    # df = df.iloc[:6]
                    # 按family_num排序
                    df = df.sort_values(by='family_num', ascending=True)
                    # 创建柱状图
                    ic(df)
                    bars = axs[n].barh(df['familyID'], df['family_num'], color=colors[n])

                # bars = axs[n].bar(df['familyID'], df['family_num'],color=colors[n])
                
                    max_val = max(df['family_num']) * 1.35 # 增加10%的缓冲区
                    axs[n].set_xlim(0, max_val)
                    for bar in bars:
                        width = bar.get_width()
                        axs[n].annotate(f'{width}',
                                        xy=(width, bar.get_y() + bar.get_height() / 2),
                                        xytext=(15, 0),  # 3 points horizontal offset
                                        textcoords="offset points",
                                        ha='center', va='center', fontsize=16, weight='bold')

                    # Set subplot titles
                    axs[n].set_title(title, fontsize=title_fontsize, weight='bold', pad=15)

                    # Set spines
                    axs[n].spines['top'].set_linewidth(1.5)
                    axs[n].spines['right'].set_linewidth(1.5)
                    axs[n].spines['left'].set_linewidth(1.5)
                    axs[n].spines['bottom'].set_linewidth(1.5)

                    # Label formatting
                    axs[n].set_ylabel('Family ID', fontsize=axis_label_fontsize, weight='bold')
                    axs[n].set_xlabel('Family Degree', fontsize=axis_label_fontsize, weight='bold')
                    axs[n].tick_params(axis='both', labelsize=18)
                    axs[n].tick_params(bottom=False, top=False, left=True, right=False)
                    axs[n].set_yticklabels(df['familyID'], fontsize=18, weight='bold')

                    for label in axs[n].get_xticklabels() + axs[n].get_yticklabels():
                        label.set_weight('bold')
                    n=n+1
            # 自动调整子图布局
        plt.tight_layout()
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/core_gene_tf/homo_{data_type}_tf_num_all3.png',bbox_inches='tight')

    #--------------------------------------------gene_num----------------------------------------------
    import pandas as pd
    from sklearn.preprocessing import MinMaxScaler
    # species= ['japonica','indica_BGI', 'zeamays']#'homo'sys.argv[1]
    # exps=['CK','CD','CE']#
    for data_type in data_types:
        # if data_type == 'KEGG':
        #     species=['indica','zeamays','japonica']
        # else:
        #     species=['indica_BGI','zeamays','japonica']
        species=['homo']
        n=0
        if 'homo' in species:
            fig, axs = plt.subplots(1, 3, figsize=(15, 5))
            colors = [(77,187,213),
                (0,160,135),(230,75,53)
                ]
        elif len(exps)==3:
            fig, axs = plt.subplots(3, 3, figsize=(15, 12))
            colors = [(77,187,213), (77,187,213), (77,187,213),
                (0,160,135),(0,160,135),(0,160,135),
                (230,75,53),(230,75,53),(230,75,53),
                ]
        else:
            fig, axs = plt.subplots(3, 2, figsize=(12, 12))
            # colors = ["#FFA500", "#FFA500","#FFA500", 
            #            "#800080", "#800080","#800080",
            #             "#008000","#008000","#008000",
            #           ]
            colors = [(77,187,213), (77,187,213),
                    (0,160,135),(0,160,135),
                    (230,75,53),(230,75,53),
                    ]
        colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
        colors= colors_normalized
        gos=['all']#,'transpiration','nitrogen','photosynthesis','water_deprivation'
        #title = [f'{exp}-all',f'{exp}-transpiration',f'{exp}-nitrogen',f'{exp}-photosynthesis',f'{exp}-water_deprivation']


        for go in gos:
            n = 0
            #fig, axs = plt.subplots(3, 2, figsize=(12, 12))
            ax = axs.flatten()
            for speci in species:
                for exp in exps:
                    csv_files = [
                        f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{exp}_c_degree_{go}.csv',
                        # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_c_closeness_{go}.csv',
                        # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_pagerank_{go}.csv',
                        # f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/{speci}/{exp}_between_{go}.csv'
                    ]

                    dataframes = []
                    for file in csv_files:
                        df = pd.read_csv(file)
                        df = df.iloc[:, [1]]  # Adjust based on actual data
                        dataframes.append(df)

                    base_df = pd.read_csv(csv_files[0]).iloc[:, [0]]
                    combined_df = pd.concat([base_df] + dataframes, axis=1)
                    scaler = MinMaxScaler(feature_range=(0, 1))
                    combined_df.iloc[:, 1:] = scaler.fit_transform(combined_df.iloc[:, 1:])
                    combined_df['avg'] = combined_df.iloc[:, 1]#[:, 1:].mean(axis=1)    只取degree部分

                    title = f'{exp}-{speci}'
                    df = combined_df.sort_values(by='avg', ascending=False)
                    if speci=='homo':
                        df = df.iloc[:6]
                    else:
                        df = df.iloc[:10]
                    df = df.sort_values(by='avg', ascending=True)
                    # Horizontal bar plot
                    bars = ax[n].barh(df.iloc[:,0], df['avg'], color=colors[n])

                    # Set axes limits
                    # if speci=='homo':
                    #     pass
                    # elif len(exps)==3:
                    #     ax[0].set_xlim(0.8,1.1)
                    #     ax[1].set_xlim(0.8,1.1)
                    #     ax[2].set_xlim(0.8,1.1)
                    #     ax[3].set_xlim(0.8,1.1)
                    #     ax[4].set_xlim(0.8,1.1)
                    #     ax[5].set_xlim(0.8,1.1)
                    #     ax[6].set_xlim(0.5,1.1)
                    #     ax[7].set_xlim(0.5,1.1)
                    #     ax[8].set_xlim(0.5,1.1)
                    # else:
                    #     ax[0].set_xlim(0.8,1.1)
                    #     ax[1].set_xlim(0.8,1.1)
                    #     ax[2].set_xlim(0.8,1.1)
                    #     ax[3].set_xlim(0.8,1.1)
                    #     ax[4].set_xlim(0.5,1.1)
                    #     ax[5].set_xlim(0.5,1.1)

                    ax[n].set_title(title, fontsize=25, weight='bold')
                    ax[n].spines['bottom'].set_linewidth(1.5)
                    ax[n].spines['right'].set_linewidth(1.5)
                    ax[n].spines['left'].set_linewidth(1.5)
                    ax[n].spines['top'].set_linewidth(1.5)
                    ax[n].tick_params(bottom=False, top=False, left=True, right=False)

                    # Label adjustments for horizontal layout
                    ax[n].set_ylabel('Gene ID', fontsize=18, weight='bold')
                    ax[n].set_xlabel('Gene Degree', fontsize=18, weight='bold')
                    ax[n].set_yticklabels(df.iloc[:,0], rotation=0, ha='right', fontsize=12, weight='bold')

                    ax[n].tick_params(axis='both', labelsize=18)
                    for label in ax[n].get_yticklabels() + ax[n].get_xticklabels():
                        label.set_weight('bold')
                    n += 1

            plt.tight_layout()
            plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/core_gene_tf/homo_{data_type}_{go}_gene_degree_all3.png', bbox_inches='tight')

conditions=['CK','CD','CE']#
data_types=['TF','PPI','KEGG']
#for data_type in data_types:

plot_core_gene(conditions,data_types)
#----------------------------------------------GRN聚类系数与平均密度------------------------------------------
# species= ['japonica','indica', 'zeamays']#sys.argv[1]
# exps=['CK','CD','CE']
# n=0
# # colors =['#556B2F','#556B2F','#556B2F','#556B2F','#556B2F',
# #         '#800000', '#800000','#800000','#800000','#800000', 
# #         '#191970','#191970','#191970','#191970','#191970']
# colors =['#800000', '#800000','#800000',
#          '#556B2F','#556B2F','#556B2F',
#         '#191970','#191970','#191970']
# gos=['all','transpiration','nitrogen','photosynthesis','water_deprivation']
# #title = [f'{exp}-all',f'{exp}-transpiration',f'{exp}-nitrogen',f'{exp}-photosynthesis',f'{exp}-water_deprivation']
# species=['indica','japonica','zeamays']
# treats=['CK','CD','CE']

# for go in gos:
#     n=0
#     fig, axs = plt.subplots(2, 3, figsize=(15, 12))
#     ax=axs.flatten()
#     all_dfs=[]
#     for speci in species:
#         dfs=[]
#         for exp in exps:
#             csv_files = f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_grn_feature_{go}.csv'
#             title = f'{exp}-{speci}'
#             title_fontsize = 20
#             axis_label_fontsize = 14
#             df =pd.read_csv(csv_files)
#             dfs.append(df['0'])
#         all_dfs.append(dfs)
#     all_dfs = [item for sublist in all_dfs for item in sublist]
#     all_dfs = np.array([s.values.tolist() for s in all_dfs])
#     cluster = all_dfs[:,0]
#     density = all_dfs[:,1]
#     all_dfs=[]
#     for speci in species:
#         dfs=[]
#         for exp in exps:
#             csv_files = f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_grn_feature_add_{go}.csv'
#             title = f'{exp}-{speci}'
#             title_fontsize = 20
#             axis_label_fontsize = 14
#             df =pd.read_csv(csv_files)
#             dfs.append(df['0'])
#         all_dfs.append(dfs)
#     # 创建柱状图
#     all_dfs = [item for sublist in all_dfs for item in sublist]
#     all_dfs = np.array([s.values.tolist() for s in all_dfs])
#     ic(all_dfs)
#     nodes = all_dfs[:,0]
#     edges = all_dfs[:,1]
#     transitivity = all_dfs[:,2]
#     reciprocity = all_dfs[:,3]
#     gene_weights=[]
#     for speci in species:
#         for exp in exps:
#             gene_weight_data = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{speci}/{exp}_gene_weights_all.csv')
#             ic(gene_weight_data['pred_weights1'].sum())
#             gene_weight=gene_weight_data['pred_weights1'].sum()
#             gene_weights.append(gene_weight)

#     # 使用列表推导式提取奇数索引的元素
#     groups = [all_dfs[i:i+3] for i in range(0, len(cluster), 3)]
#     positions = range(len(cluster))
#     all_data=[transitivity,cluster,gene_weights,nodes,edges,density]
#     titles=['transitivity','cluster','gene-weights','nodes','edges','density']
#     ic(cluster)
#     for data in all_data:

#         ax[n].tick_params(axis='x', labelsize=22)
#         ax[n].tick_params(axis='x', labelsize=22)
#         ax[n].bar(positions, data, color=colors)  # 你可以选择不同的颜色z
#         ax[n].set_xticks(positions, ['indica CK','indica CD','indica CE','japonica CK','japonica CD','japonica CE','zeamays CK','zeamays CD','zeamays CE'],rotation=90)  # 创建标签列表匹配自定义位置
#         # 设置图表的标题和坐标轴标签
#         ax[n].set_title(f'{titles[n]}',fontsize=30,weight='bold')
#         ax[n].set_xlabel('Group',fontsize=25,weight='bold')
#         ax[n].set_ylabel(f'{titles[n]}',fontsize=25,weight='bold')
#         ax[n].spines['top'].set_linewidth(1.5)
#         ax[n].spines['bottom'].set_linewidth(1.5)
#         ax[n].spines['right'].set_linewidth(1.5)
#         ax[n].spines['left'].set_linewidth(1.5)
#         # current_labels = [item.get_text() for item in ax[n].get_xticklabels()]
#         # ax[n].set_xticklabels(current_labels,df.iloc[:,0], rotation=90, ha='right', fontsize=12,weight='bold')
#         # current_labels = [item.get_text() for item in ax[n].get_yticklabels()]
#         # ax[n].set_yticklabels(current_labels,rotation=90, ha='right', fontsize=12,weight='bold')
#         #ax[n].tick_params(axis='both', labelsize=14)
#         ax[n].tick_params(axis='x', labelrotation=90, labelright=True, labelsize=22)
#         ax[n].tick_params(axis='y',labelsize=22)
#         for label in ax[n].get_xticklabels() + ax[n].get_yticklabels():
#             label.set_weight('bold')
#         current_axes = ax[1]  # 获取当前轴对象
#         spines_to_adjust = ['top', 'bottom', 'left', 'right']  # 定义要调整的边框列表
#         for spine in spines_to_adjust:
#             current_axes.spines[spine].set_linewidth(2)
#         n=n+1

#     # 显示图表
#     plt.tight_layout(pad=1.0)
#     plt.savefig(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{go}_RGN_feature.png',bbox_inches='tight')

# #=================================================family weights============================================
# fig, axs = plt.subplots(5, 1, figsize=(15, 12))
# csv_files = [f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_data.txt',f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_transpiration.csv', 
#             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_nitrogen.csv', 
#             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_photosynthesis.csv',
#              f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights_water_deprivation.csv']
# title = [f'{exp}-transpiration',f'{exp}-nitrogen',f'{exp}-photosynthesis',f'{exp}-water_deprivation']
# # 遍历CSV文件和子图轴
# # 定义标题和轴标签的字体大小
# title_fontsize = 20
# axis_label_fontsize = 14

# i=0
# for ax, csv_file in zip(axs.flatten(), csv_files):
#     # 加载CSV文件
#     df = pd.read_csv(csv_file)
#     df = df.sort_values(by='pred_weights1', ascending=False)
#     # 删除第一行
#     df = df.iloc[1:]

#     # 按family_num排序
    
#     df = df.iloc[:50]
#     # 创建柱状图
#     bars=ax.bar(df['familyID'], df['pred_weights1'])
#     for bar in bars:
#         height = bar.get_height()
#         ax.annotate(f'{height:.0f}',
#                 xy=(bar.get_x() + bar.get_width() / 2, height),
#                 xytext=(0, 3),  # 3点垂直偏移
#                 textcoords="offset points",
#                 ha='center', va='bottom', fontsize=12, weight='bold')
#     # 设置子图标题
#     ax.set_title(f'{title[i]}',fontsize=title_fontsize,weight='bold',)
#     ax.spines['bottom'].set_linewidth(1.5)  # 设置边框线宽为2.0
#     ax.spines['right'].set_linewidth(1.5)
#     ax.spines['left'].set_linewidth(1.5)
#     ax.spines['top'].set_linewidth(1.5)
#     ax.tick_params(bottom=False, top=False, left=True, right=False)
#     # 设置坐标轴标签
#     ax.set_xlabel('Family ID', fontsize=axis_label_fontsize,weight='bold')
#     ax.set_ylabel('Family Weights', fontsize=axis_label_fontsize,weight='bold')
#     #ax.set_xticks(range(len(df['familyID'])))
    
#     ax.set_xticklabels(df['familyID'], rotation=60, ha='right', fontsize=12,weight='bold')
#     #ax.set_yticklabels( fontsize=12,weight='bold')
#     ax.tick_params(axis='both', labelsize=14)
#     for label in ax.get_xticklabels() + ax.get_yticklabels():
#         label.set_weight('bold')
#     #ax.tick_params(axis='y', labelsize=12, labelweight='bold')
#     #ax.set_yticks(ticks=[],fontproperties='Times New Roman', size=14)
#     #ax.set_xticks(ticks=[],fontproperties='Times New Roman', size=14)
#     i=i+1
# # 自动调整子图布局
# plt.tight_layout()
# plt.savefig(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_family_weights.png',bbox_inches='tight')


# #=============================================gene weights========================================


# fig, axs = plt.subplots(5, 1, figsize=(15, 12))
# csv_files = [f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights_data.txt',f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights_transpiration.csv', 
#             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights_nitrogen.csv', 
#             f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights_photosynthesis.csv',
#              f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights_water_deprivation.csv']
# title = [f'{exp}-transpiration',f'{exp}-nitrogen',f'{exp}-photosynthesis',f'{exp}-water_deprivation']
# # 遍历CSV文件和子图轴
# # 定义标题和轴标签的字体大小
# title_fontsize = 20
# axis_label_fontsize = 14

# i=0
# for ax, csv_file in zip(axs.flatten(), csv_files):
#     # 加载CSV文件
#     df = pd.read_csv(csv_file)

#     # 删除第一行
#     #df = df.iloc[1:]

#     # 按family_num排序
#     df = df.sort_values(by='pred_weights1', ascending=False)
#     df = df.iloc[:50]
#     # 创建柱状图
#     bars=ax.bar(df['GeneID'], df['pred_weights1'])
#     for bar in bars:
#         height = bar.get_height()
#         ax.annotate(f'{height:.0f}',
#                         xy=(bar.get_x() + bar.get_width() / 2, height),
#                         xytext=(0, 3),  # 3点垂直偏移
#                         textcoords="offset points",
#                         ha='center', va='bottom', fontsize=12, weight='bold')
#     # 设置子图标题
#     ax.set_title(f'{title[i]}',fontsize=title_fontsize,weight='bold',)
#     ax.spines['bottom'].set_linewidth(1.5)  # 设置边框线宽为2.0
#     ax.spines['right'].set_linewidth(1.5)
#     ax.spines['left'].set_linewidth(1.5)
#     ax.spines['top'].set_linewidth(1.5)
#     ax.tick_params(bottom=False, top=False, left=True, right=False)
#     # 设置坐标轴标签
#     ax.set_xlabel('Gene ID', fontsize=axis_label_fontsize,weight='bold')
#     ax.set_ylabel('Gene Weights', fontsize=axis_label_fontsize,weight='bold')
#     #ax.set_xticks(range(len(df['familyID'])))
    
#     ax.set_xticklabels(df['GeneID'], rotation=60, ha='right', fontsize=12,weight='bold')
#     #ax.set_yticklabels( fontsize=12,weight='bold')
#     ax.tick_params(axis='both', labelsize=14)
#     for label in ax.get_xticklabels() + ax.get_yticklabels():
#         label.set_weight('bold')
#     #ax.tick_params(axis='y', labelsize=12, labelweight='bold')
#     #ax.set_yticks(ticks=[],fontproperties='Times New Roman', size=14)
#     #ax.set_xticks(ticks=[],fontproperties='Times New Roman', size=14)
#     i=i+1
# # 自动调整子图布局
# plt.tight_layout()
# plt.savefig(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/{exp}_gene_weights.png',bbox_inches='tight')