import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import igraph as ig
from icecream import install,ic
install()
config = {
    "font.family":'Times New Roman',
    "font.size": 18,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
}
plt.rcParams.update(config)


# data_types=['TF','PPI','KEGG']
# for data_type in data_types:
#     if data_type== 'KEGG':
#         species=['indica','zeamays','japonica']
#     else:
#         species=['indica_BGI','zeamays','japonica']

# #species=['indica_BGI','japonica','zeamays']#'homo'
#     titles=['indica','japonica','zeamays']
#     treats=['CK','CD']#,'CE'
#     all_weights=[]
#     #layers=['x2']#,'x3''x0','x1',
#     organs=['leaf','root']
#     #axs= plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2, fig=plt.gcf())
#     # colors = [(77,187,213), (230,75,53),(0,160,135)]
#     # colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
#     # colors=colors_normalized
#     colors = ['black','red']

#     # colors = [(0.7, 0.1, 0.3),  # 明亮红
#     #           (0.1, 0.5, 0.2),  # 明亮绿
#     #           (0.1, 0.1, 0.7)]  # 明亮蓝

#     #for layer in layers:


    # gos=['nitrogen','photosynthesis','stoma','transpiration','water_deprivation']#
    # fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    # axs=axs.flatten()
    
    # for go in gos:
    #     i = 0
    #     for n, speci in enumerate(species):
    #         avgds = []
    #         avg_clusters = []
            
    #         for treat in treats:
    #             file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_{go}_rows.csv'
    #             data = pd.read_csv(file_path)
    #             G = nx.DiGraph()
    #             G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
                
    #             d = dict(nx.degree(G, weight='pred_weights2'))  # 平均度
    #             avg_d = sum(d.values()) / len(G.nodes)
    #             ic(f'{speci}_{treat}_{go}', avg_d)
    #             ic(f'{speci}_{treat}_{go}', len(G.nodes))
                
    #             G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
    #             avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
    #             avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
                
    #             avgds.append(avg_d)
    #             avg_clusters.append(avg_cluster)
            
    #         labels = ['CK', 'CD']
    #         bar_width = 0.4
    #         x = np.arange(len(labels))
            
    #         # 绘制平均度的柱状图
    #         axs[i].bar(x - bar_width / 2, avgds, bar_width, label='Average Degree', color=['black', 'red'])
    #         for j, val in enumerate(avgds):
    #             axs[i].text(j - bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')
            
    #         # 绘制平均聚类系数的柱状图
    #         axs[i].bar(x + bar_width / 2, avg_clusters, bar_width, label='Cluster Coefficient', color=['grey', 'blue'])
    #         for j, val in enumerate(avg_clusters):
    #             axs[i].text(j + bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')
            
    #         axs[i].set_xticks(x)
    #         axs[i].set_xticklabels(labels, fontsize=22, rotation=0, weight='bold')
            
    #         y_min = min(min(avgds), min(avg_clusters)) - (max(max(avgds), max(avg_clusters)) - min(min(avgds), min(avg_clusters))) * 0.1
    #         y_max = max(max(avgds), max(avg_clusters)) + (max(max(avgds), max(avg_clusters)) - min(min(avgds), min(avg_clusters))) * 0.1
    #         axs[i].set_ylim(y_min, y_max)
            
    #         axs[i].set_title(f'{titles[i]}', fontsize=30, weight='bold', pad=15)
    #         axs[i].set_xlabel('Group', fontsize=25, weight='bold')
    #         axs[i].set_ylabel('Values', fontsize=25, weight='bold')
    #         axs[i].spines['top'].set_linewidth(1.5)
    #         axs[i].spines['bottom'].set_linewidth(1.5)
    #         axs[i].spines['right'].set_linewidth(1.5)
    #         axs[i].spines['left'].set_linewidth(1.5)
    #         axs[i].legend(fontsize=15)
            
    #         i += 1

    #     plt.tight_layout()
    #     plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_{go}_dgree2_sink.png', bbox_inches='tight')
    #     plt.clf()

    # fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    # axs=axs.flatten()
    # i=0
    # for n, speci in enumerate(species):
    #     avgds = []
    #     avg_clusters = []

    #     for treat in treats:
    #         file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_{go}_rows.csv'
    #         data = pd.read_csv(file_path)
    #         G = nx.DiGraph()
    #         G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
            
    #         d = dict(nx.degree(G, weight='pred_weights2'))  # 平均度
    #         avg_d = sum(d.values()) / len(G.nodes)
    #         ic(f'{speci}_{treat}_{go}', avg_d)
    #         ic(f'{speci}_{treat}_{go}', len(G.nodes))
            
    #         G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
    #         avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
    #         avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
            
    #         avgds.append(avg_d)
    #         avg_clusters.append(avg_cluster)
        
    #     labels = ['CK', 'CD']
    #     bar_width = 0.4
    #     x = np.arange(len(labels))
        
    #     # 绘制平均度的柱状图
    #     axs[i].bar(x - bar_width / 2, avgds, bar_width, label='Average Degree', color=['black', 'red'])
    #     for j, val in enumerate(avgds):
    #         axs[i].text(j - bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')
        
    #     # 绘制平均聚类系数的柱状图
    #     axs[i].bar(x + bar_width / 2, avg_clusters, bar_width, label='Cluster Coefficient', color=['grey', 'blue'])
    #     for j, val in enumerate(avg_clusters):
    #         axs[i].text(j + bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')
        
    #     axs[i].set_xticks(x)
    #     axs[i].set_xticklabels(labels, fontsize=22, rotation=0, weight='bold')
        
    #     y_min = min(min(avgds), min(avg_clusters)) - (max(max(avgds), max(avg_clusters)) - min(min(avgds), min(avg_clusters))) * 0.1
    #     y_max = max(max(avgds), max(avg_clusters)) + (max(max(avgds), max(avg_clusters)) - min(min(avgds), min(avg_clusters))) * 0.1
    #     axs[i].set_ylim(y_min, y_max)
        
    #     axs[i].set_title(f'{titles[i]}', fontsize=30, weight='bold', pad=15)
    #     axs[i].set_xlabel('Group', fontsize=25, weight='bold')
    #     axs[i].set_ylabel('Values', fontsize=25, weight='bold')
    #     axs[i].spines['top'].set_linewidth(1.5)
    #     axs[i].spines['bottom'].set_linewidth(1.5)
    #     axs[i].spines['right'].set_linewidth(1.5)
    #     axs[i].spines['left'].set_linewidth(1.5)
    #     axs[i].legend(fontsize=15)
        
    #     i += 1

    # plt.tight_layout()
    # plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_allGRN_dgree2_sink.png', bbox_inches='tight')
    # plt.clf()

#================================================加上复水的结果===================================================================
# data_types=['TF','PPI','KEGG']
# for data_type in data_types:
#     if data_type== 'KEGG':
#         species=['indica','zeamays','japonica']
#     else:
#         species=['indica_BGI','zeamays','japonica']

# #species=['indica_BGI','japonica','zeamays']#'homo'
#     titles=['indica','japonica','zeamays']
#     treats=['CK','CD','CE']#
#     all_weights=[]
#     #layers=['x2']#,'x3''x0','x1',
#     organs=['leaf','root']
#     #axs= plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2, fig=plt.gcf())
#     # colors = [(77,187,213), (230,75,53),(0,160,135)]
#     # colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
#     # colors=colors_normalized
#     colors = ['black','red','grey']

#     # colors = [(0.7, 0.1, 0.3),  # 明亮红
#     #           (0.1, 0.5, 0.2),  # 明亮绿
#     #           (0.1, 0.1, 0.7)]  # 明亮蓝

#     #for layer in layers:


#     gos=['nitrogen','photosynthesis','stoma','transpiration','water_deprivation']#
#     #fig, axs = plt.subplots(1, 3, figsize=(15, 6))
#     #axs=axs.flatten()
#     for go in gos:
#         colors = ['black', 'red', 'grey']
#         labels = ['CK', 'CD', 'CE']
#         fig, ax = plt.subplots(figsize=(8, 5))
#         labels_plot = ['indica', 'japonica', 'zeamays']
#         bar_width = 0.4
#         gap_width = 0.8
#         tick_positions = []

#         # 循环处理每个物种的数据
#         for i, speci in enumerate(species):
#             avgds = []
#             avg_clusters = []
#             for treat in treats:
#                 file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_{go}_rows.csv'
#                 data = pd.read_csv(file_path)
#                 # if data.empty:
#                 #     data =[]
#                 G = nx.DiGraph()
#                 G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
                
#                 d = dict(nx.degree(G, weight='pred_weights2'))  # 平均度
#                 avg_d = sum(d.values()) / len(G.nodes)
#                 ic(f'{speci}_{treat}_{go}', avg_d)
#                 ic(f'{speci}_{treat}_{go}', len(G.nodes))
                
#                 G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
#                 avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
#                 avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
                
#                 avgds.append(avg_d)
#                 avg_clusters.append(avg_cluster)

#             x = np.arange(len(labels)) + i * (len(labels) + gap_width)  # 计算每组柱的位置
#             print(len(avgds))
#             ax.bar(x, avgds, width=bar_width, color=colors, label=labels_plot[i])  # 绘制柱状图

#             tick_positions.append(x.mean())  # 计算每组柱的中间位置

#             for j, (treat, color) in enumerate(zip(treats, colors)):
#                 ax.bar(x[j], avgds[j], width=bar_width, color=color, label=treat if i == 0 else "")  # 绘制柱状图
#                 #ax.text(x[j], avgds[j], f'{avgds[j]:.3f}', ha='center', va='bottom', fontsize=16, weight='bold')
#             for j, val in enumerate(avgds):
#                 ax.text(x[j], val, f'{val:.3f}', ha='center', va='bottom', fontsize=16, weight='bold')

#         # 设置x轴标签和标题
#         ax.set_xticks(tick_positions)
#         ax.set_xticklabels(labels_plot, fontsize=22, rotation=0, weight='bold')
#         ax.set_title(f'{go} combined degree', fontsize=30, weight='bold', pad=15)
#         ax.set_ylabel('Average Degree', fontsize=25, weight='bold')
#         ax.spines['top'].set_linewidth(1.5)
#         ax.spines['bottom'].set_linewidth(1.5)
#         ax.spines['right'].set_linewidth(1.5)
#         ax.spines['left'].set_linewidth(1.5)

#         # 设置图例
#         handles, _ = ax.get_legend_handles_labels()
#         by_label = dict(zip(labels, handles))
#         ax.legend(by_label.values(), by_label.keys(),  fontsize=16, title_fontsize='17')
#         #ax.legend()

#         plt.tight_layout()
#         plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_{go}_combined_degree_sink.png', bbox_inches='tight')
#         plt.clf()
import matplotlib.patches as mpatches
data_types = ['TF', 'PPI', 'KEGG']
for data_type in data_types:
    if data_type == 'KEGG':
        species = ['indica', 'zeamays', 'japonica']
    else:
        species = ['indica_BGI', 'zeamays', 'japonica']
        
    labels = ['CK', 'CD', 'CE']
    titles = ['indica', 'japonica', 'zeamays']
    treats = ['CK', 'CD', 'CE']
    colors = ['black', 'red', 'grey']
    
    fig, ax = plt.subplots(figsize=(9, 5))
    labels_plot = ['indica', 'japonica', 'zeamays']
    
    bar_width = 0.2  # 单个柱子的宽度
    gap_width = 0.4  # 每个物种组之间的间隔宽度
    tick_positions = []
    
    for n, speci in enumerate(species):
        avgds = []
        avg_clusters = []
        
        for treat in treats:
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')
            G = nx.DiGraph()
            G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
            
            d = dict(nx.degree(G, weight='pred_weights2'))  # 平均度
            avg_d = sum(d.values()) / len(G.nodes)
            G_igraph = ig.Graph.from_networkx(G)
            avg_cluster = sum(G_igraph.transitivity_local_undirected(mode="zero")) / len(G_igraph.vs)
            
            avgds.append(avg_d)
            avg_clusters.append(avg_cluster)
        
        x = np.arange(len(treats)) * (bar_width + 0.1) + n * (bar_width * len(treats) + gap_width)  # 计算每组柱子的位置
        ax.bar(x, avgds, width=bar_width, color=colors[:len(treats)], label=labels_plot[n])
        
        tick_positions.append(x.mean())  # 计算每组柱的中间位置

        for j, val in enumerate(avgds):
            ax.text(x[j], val, f'{val:.3f}', ha='center', va='bottom', fontsize=16, weight='bold')

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(labels_plot, fontsize=22, rotation=0, weight='bold')
    ax.set_title(f'Combined degree', fontsize=30, weight='bold', pad=15)
    ax.set_ylabel('Average Degree', fontsize=25, weight='bold')
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)

    # 设置图例
    legend_colors = ['black', 'red', 'grey']  # 这里可以调整颜色以适应不同的需求
    legend_labels = ['CK', 'CD', 'CE']

    # 创建图例标签
    legend_handles = [mpatches.Patch(color=color, label=label) for color, label in zip(legend_colors, legend_labels)]
    ax.legend(handles=legend_handles, fontsize=16, title_fontsize='17')

    plt.tight_layout()
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_all_combined_degree_sink.png', bbox_inches='tight')
    plt.clf()

        #colors = ['black', 'red', 'grey']

        # fig, ax = plt.subplots(nrows=2,ncols=1,figsize=(12, 6))
        # for i, speci in enumerate(species):
        #     avgds = []
        #     avg_clusters = []

        #     for treat in treats:
        #         file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
        #         data = pd.read_csv(file_path, sep='\t')
        #         G = nx.DiGraph()
        #         G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
                
        #         d = dict(nx.degree(G, weight='pred_weights2'))
        #         avg_d = sum(d.values()) / len(G.nodes)
        #         ic(f'{speci}_{treat}_{go}', avg_d)
        #         ic(f'{speci}_{treat}_{go}', len(G.nodes))

        #         G_igraph = ig.Graph.from_networkx(G)
        #         avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")
        #         avg_cluster = sum(avg_cluster) / len(avg_cluster)

        #         avgds.append(avg_d)
        #         avg_clusters.append(avg_cluster)

        #     bar_width = 0.4
        #     x = np.arange(len(treats))

        #     # 绘制平均度的柱状图
        #     ax.bar(x - bar_width / 2, avgds, bar_width, label='Average Degree', color=colors)
        #     for j, val in enumerate(avgds):
        #         ax.text(j - bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')

        #     # 绘制平均聚类系数的柱状图
        #     ax.bar(x + bar_width / 2, avg_clusters, bar_width, label='Cluster Coefficient', color=colors)
        #     for j, val in enumerate(avg_clusters):
        #         ax.text(j + bar_width / 2, val, f'{val:.3f}', ha='center', va='bottom', fontsize=12, weight='bold')

        #     ax.set_xticks(x)
        #     ax.set_xticklabels(treats, fontsize=22, rotation=0, weight='bold')

        #     ax.set_title(f'{titles[i]}', fontsize=30, weight='bold', pad=15)
        #     ax.set_xlabel('Group', fontsize=25, weight='bold')
        #     ax.set_ylabel('Values', fontsize=25, weight='bold')
        #     ax.spines['top'].set_linewidth(1.5)
        #     ax.spines['bottom'].set_linewidth(1.5)
        #     ax.spines['right'].set_linewidth(1.5)
        #     ax.spines['left'].set_linewidth(1.5)
        #     ax.legend(fontsize=15)

        # plt.tight_layout()
        # plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_allGRN_dgree3_sink.png', bbox_inches='tight')
        # plt.clf()






