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


data_types=['TF','PPI','KEGG']#
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']

#species=['indica_BGI','japonica','zeamays']#'homo'
    titles=['indica','japonica','zeamays']
    treats=['CK','CD']#,'CE'
    all_weights=[]
    #layers=['x2']#,'x3''x0','x1',
    organs=['leaf','root']
    #axs= plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2, fig=plt.gcf())
    # colors = [(77,187,213), (230,75,53),(0,160,135)]
    # colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
    # colors=colors_normalized
    colors = ['black','red']

    # colors = [(0.7, 0.1, 0.3),  # 明亮红
    #           (0.1, 0.5, 0.2),  # 明亮绿
    #           (0.1, 0.1, 0.7)]  # 明亮蓝

    #for layer in layers:


    gos=['nitrogen','photosynthesis','stoma','transpiration','water_deprivation']#
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    axs=axs.flatten()

    for go in gos:
        i=0
        fig, axs = plt.subplots(1, 3, figsize=(15, 6)) 
        axs=axs.flatten()
        for n,speci in enumerate(species):
            avgds=[]
            d_hiss=[]
            avg_clusters=[]
            #for organ in organs:
            for treat in treats:
                file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_{go}_rows.csv'
                data = pd.read_csv(file_path)
                G = nx.DiGraph()
                #G=nx.from_pandas_edgelist(data, 'fromNode', 'toNode',edge_attr='weight',create_using=nx.DiGraph())
                G=nx.from_pandas_edgelist(data, 'TF', 'Target',edge_attr='pred_weights2',create_using=nx.DiGraph())
                
                d = dict(nx.degree(G,weight='pred_weights2'))  #平均度
                avg_d = sum(d.values())/len(G.nodes)
                ic(f'{speci}_{treat}_{go}',avg_d)
                ic(f'{speci}_{treat}_{go}',len(G.nodes))
                #d_his=nx.degree_histogram(G) #度分布
                #ic(f'{speci}_{treat}',d_his)
                
                G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
                avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
                avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
                # if speci=='indica' and organ=='leaf' and treat=='CD':
                #     avg_cluster=0.709452335236
                #avg_cluster=nx.average_clustering(G,weight='weight')
                #ic(f'{speci}_{treat}',avg_cluster)
                avgds.append(avg_d)
                #d_hiss.append(d_his)
                avg_clusters.append(avg_cluster)
            #labels=['leaf CK','leaf CD','root CK','root CD']
            labels=['CK','CD']
            axs[i].bar(range(len(avgds)), avgds, align='center',color=['black','red','black','red']) # 假设d是一个字典，转换为值列表进行绘制
            axs[i].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i].set_xticklabels(labels,fontsize=22,rotation=0,weight='bold')
            # max_val = max(avgds) * 1.1  # 增加10%的缓冲区
            # axs[i].set_ylim(0, max_val)
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i].set_ylim(y_min, y_max)
            for j, val in enumerate(avgds):
                    axs[i].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold') 


            axs[i].set_title(f'{titles[i]} Degree',fontsize=30,weight='bold',pad=15)
            axs[i].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i].set_ylabel('Average Degree',fontsize=25,weight='bold')
            axs[i].spines['top'].set_linewidth(1.5)
            axs[i].spines['bottom'].set_linewidth(1.5)
            axs[i].spines['right'].set_linewidth(1.5)
            axs[i].spines['left'].set_linewidth(1.5)
            # if go=='nitrogen':             
            #     axs[0].set_ylim(2.5,2.8)
            #     axs[1].set_ylim(8,8.8)
            #     axs[2].set_ylim(2.9,3.1)
            # elif go=='photosynthesis':
            #     axs[0].set_ylim(3.2,3.7)
            #     axs[1].set_ylim(2.6,3.2)
            #     axs[2].set_ylim(3.2,3.7)             
            # elif go=='stoma':
            #     axs[0].set_ylim(2.0,2.3)
            #     axs[1].set_ylim(2.9,3.0)
            #     axs[2].set_ylim(2,2.3)         
            # elif go=='transpiration':
            #     axs[0].set_ylim(2.0,2.2)
            #     axs[1].set_ylim(2.1,2.3)
            #     axs[2].set_ylim(1.5,1.8)    
            # elif go=='water_deprivation':
            #     axs[0].set_ylim(3.0,3.2)
            #     axs[1].set_ylim(3.4,3.5)
            #     axs[2].set_ylim(1.9,2.0)  
            i=i+1
        plt.tight_layout()
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_{go}_dgree2.png',bbox_inches='tight')
        plt.clf()


    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs=axs.flatten()
    i=0
    for n,speci in enumerate(species):
        avgds=[]
        d_hiss=[]
        avg_clusters=[]
        #for organ in organs:
        for treat in treats:
            #file_path = f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_WGCNA/{speci}/out/{treat}/network_edge.txt'
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')
            G = nx.DiGraph()
            #G=nx.from_pandas_edgelist(data, 'fromNode', 'toNode',edge_attr='weight',create_using=nx.DiGraph())
            G=nx.from_pandas_edgelist(data, 'TF', 'Target',edge_attr='pred_weights2',create_using=nx.DiGraph())
            
            d = dict(nx.degree(G,weight='pred_weights2'))  #平均度
            avg_d = sum(d.values())/len(G.nodes)
            ic(f'{speci}_{treat}_{go}',avg_d)
            ic(f'{speci}_{treat}_{go}',len(G.nodes))
            #d_his=nx.degree_histogram(G) #度分布
            #ic(f'{speci}_{treat}',d_his)
            
            G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
            avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
            avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
            avgds.append(avg_d)
            #d_hiss.append(d_his)
            avg_clusters.append(avg_cluster)
        #labels=['leaf CK','leaf CD','root CK','root CD']
        labels=['CK','CD']
        axs[i].bar(range(len(avgds)), avgds, align='center',color=['black','red','black','red']) # 假设d是一个字典，转换为值列表进行绘制
        y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
        y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
        axs[i].set_ylim(y_min, y_max)
        axs[i].set_xticks(range(len(labels)))
        # 设置x轴的刻度标签
        axs[i].set_xticklabels(labels,fontsize=22,rotation=0,weight='bold')
        # max_val = max(avgds) * 1.1  # 增加10%的缓冲区
        # axs[i].set_ylim(0, max_val)
        for j, val in enumerate(avgds):
                axs[i].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold') 

        axs[i].set_title(f'{titles[i]} Degree',fontsize=30,weight='bold',pad=15)
        axs[i].set_xlabel('Group',fontsize=25,weight='bold')
        axs[i].set_ylabel('Average Degree',fontsize=25,weight='bold')
        axs[i].spines['top'].set_linewidth(1.5)
        axs[i].spines['bottom'].set_linewidth(1.5)
        axs[i].spines['right'].set_linewidth(1.5)
        axs[i].spines['left'].set_linewidth(1.5)
        #第三个子图：平均聚类系数的柱状图
        #注意，avg_cluster应该只有一个值，如果它是一个列表且您想显示每个treat的结果，则需要调整
        if speci=='homo':
            axs[i+1].bar(range(len(avg_clusters)), avg_clusters,color=['black','red','black','red'])
            # max_val = max(avg_clusters) * 1.1  # 增加10%的缓冲区
            # axs[i+1].set_ylim(0, max_val)
            avgds = avg_clusters
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i+1].set_ylim(y_min, y_max)
            for j, val in enumerate(avg_clusters):
                axs[i+1].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold')

            axs[i+1].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i+1].set_xticklabels(labels,fontsize=22,weight='bold',rotation=0)
            axs[i+1].set_title(f'{titles[i]} Cluster Coef',fontsize=30,weight='bold',pad=15)
            axs[i+1].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i+1].set_ylabel('Coefficient',fontsize=25,weight='bold')
            axs[i+1].spines['top'].set_linewidth(1.5)
            axs[i+1].spines['bottom'].set_linewidth(1.5)
            axs[i+1].spines['right'].set_linewidth(1.5)
            axs[i+1].spines['left'].set_linewidth(1.5)
        else:
            axs[i+3].bar(range(len(avg_clusters)), avg_clusters,color=['black','red','black','red'])
            # max_val = max(avg_clusters) * 1.1  # 增加10%的缓冲区
            # axs[i+3].set_ylim(0, max_val)
            avgds = avg_clusters
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i+3].set_ylim(y_min, y_max)
            for j, val in enumerate(avg_clusters):
                axs[i + 3].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold')

            axs[i+3].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i+3].set_xticklabels(labels,fontsize=22,weight='bold',rotation=0)
            axs[i+3].set_title(f'{titles[i]} Cluster Coef',fontsize=30,weight='bold',pad=15)
            axs[i+3].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i+3].set_ylabel('Coefficient',fontsize=25,weight='bold')
            axs[i+3].spines['top'].set_linewidth(1.5)
            axs[i+3].spines['bottom'].set_linewidth(1.5)
            axs[i+3].spines['right'].set_linewidth(1.5)
            axs[i+3].spines['left'].set_linewidth(1.5)
            # axs[0].set_ylim(23,25)
            # axs[1].set_ylim(27,28.5)
            # axs[2].set_ylim(12,13)
            # axs[3].set_ylim(0.6,0.7)
            # axs[4].set_ylim(0.5,0.6)
            # axs[5].set_ylim(0.3,0.4)
        i=i+1

    plt.tight_layout()
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_allGRN_dgree2.png',bbox_inches='tight')


#================================================加上复水的结果===================================================================
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']

#species=['indica_BGI','japonica','zeamays']#'homo'
    titles=['indica','japonica','zeamays']
    treats=['CK','CD','CE']#
    all_weights=[]
    #layers=['x2']#,'x3''x0','x1',
    organs=['leaf','root']
    #axs= plt.subplot2grid((3, 3), (0, 0), colspan=2, rowspan=2, fig=plt.gcf())
    # colors = [(77,187,213), (230,75,53),(0,160,135)]
    # colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]
    # colors=colors_normalized
    colors = ['black','red','grey']

    # colors = [(0.7, 0.1, 0.3),  # 明亮红
    #           (0.1, 0.5, 0.2),  # 明亮绿
    #           (0.1, 0.1, 0.7)]  # 明亮蓝

    #for layer in layers:


    gos=['nitrogen','photosynthesis','stoma','transpiration','water_deprivation']#
    fig, axs = plt.subplots(1, 3, figsize=(15, 6))
    axs=axs.flatten()
    for go in gos:
        i=0
        fig, axs = plt.subplots(1, 3, figsize=(15, 6)) 
        axs=axs.flatten()
        for n,speci in enumerate(species):
            avgds=[]
            d_hiss=[]
            avg_clusters=[]
            for treat in treats:
                file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/sub_GO_module/{data_type}_{speci}_{treat}_gene_{go}_rows.csv'
                data = pd.read_csv(file_path)
                if data.empty:
                    continue
                G = nx.DiGraph()
                #G=nx.from_pandas_edgelist(data, 'fromNode', 'toNode',edge_attr='weight',create_using=nx.DiGraph())
                G=nx.from_pandas_edgelist(data, 'TF', 'Target',edge_attr='pred_weights2',create_using=nx.DiGraph())
                
                d = dict(nx.degree(G,weight='pred_weights2'))  #平均度
                avg_d = sum(d.values())/len(G.nodes)
                ic(f'{speci}_{treat}_{go}',avg_d)
                ic(f'{speci}_{treat}_{go}',len(G.nodes))
                #d_his=nx.degree_histogram(G) #度分布
                #ic(f'{speci}_{treat}',d_his)
                
                G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
                avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
                avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
                # if speci=='indica' and organ=='leaf' and treat=='CD':
                #     avg_cluster=0.709452335236
                #avg_cluster=nx.average_clustering(G,weight='weight')
                #ic(f'{speci}_{treat}',avg_cluster)
                avgds.append(avg_d)
                #d_hiss.append(d_his)
                avg_clusters.append(avg_cluster)
            #labels=['leaf CK','leaf CD','root CK','root CD']
            labels=['CK','CD','CE']
            axs[i].bar(range(len(avgds)), avgds, align='center',color=['black','red','grey','black','red','grey']) # 假设d是一个字典，转换为值列表进行绘制
            #avgds = avg_clusters
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i].set_ylim(y_min, y_max)
            axs[i].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i].set_xticklabels(labels,fontsize=22,rotation=0,weight='bold')
            # max_val = max(avgds) * 1.1  # 增加10%的缓冲区
            # axs[i].set_ylim(0, max_val)
            for j, val in enumerate(avgds):
                    axs[i].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold') 


            axs[i].set_title(f'{titles[i]} Degree',fontsize=30,weight='bold',pad=15)
            axs[i].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i].set_ylabel('Average Degree',fontsize=25,weight='bold')
            axs[i].spines['top'].set_linewidth(1.5)
            axs[i].spines['bottom'].set_linewidth(1.5)
            axs[i].spines['right'].set_linewidth(1.5)
            axs[i].spines['left'].set_linewidth(1.5)

            i=i+1
        plt.tight_layout()
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_{go}_dgree3.png',bbox_inches='tight')
        plt.clf()

    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs=axs.flatten()
    i=0
    for n,speci in enumerate(species):
        avgds=[]
        d_hiss=[]
        avg_clusters=[]
        #for organ in organs:
        for treat in treats:
            #file_path = f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_WGCNA/{speci}/out/{treat}/network_edge.txt'
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')
            G = nx.DiGraph()
            #G=nx.from_pandas_edgelist(data, 'fromNode', 'toNode',edge_attr='weight',create_using=nx.DiGraph())
            G=nx.from_pandas_edgelist(data, 'TF', 'Target',edge_attr='pred_weights2',create_using=nx.DiGraph())
            
            d = dict(nx.degree(G,weight='pred_weights2'))  #平均度
            avg_d = sum(d.values())/len(G.nodes)
            ic(f'{speci}_{treat}_{go}',avg_d)
            ic(f'{speci}_{treat}_{go}',len(G.nodes))
            #d_his=nx.degree_histogram(G) #度分布s
            #ic(f'{speci}_{treat}',d_his)
            
            G_igraph = ig.Graph.from_networkx(G)  # 将NetworkX图转换为igraph图
            avg_cluster = G_igraph.transitivity_local_undirected(mode="zero")  # 计算局部聚类系数，注意根据实际情况调整参数
            avg_cluster = sum(avg_cluster) / len(avg_cluster)  # 计算平均值
            avgds.append(avg_d)
            #d_hiss.append(d_his)
            avg_clusters.append(avg_cluster)
        #labels=['leaf CK','leaf CD','root CK','root CD']
        labels=['CK','CD','CE']
        axs[i].bar(range(len(avgds)), avgds, align='center',color=['black','red','grey','black','red','grey']) # 假设d是一个字典，转换为值列表进行绘制
        #avgds = avg_clusters
        y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
        y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
        axs[i].set_ylim(y_min, y_max)
        axs[i].set_xticks(range(len(labels)))
        # 设置x轴的刻度标签
        axs[i].set_xticklabels(labels,fontsize=22,rotation=0,weight='bold')
        # max_val = max(avgds) * 1.1  # 增加10%的缓冲区
        # axs[i].set_ylim(0, max_val)
        for j, val in enumerate(avgds):
                axs[i].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold') 

        axs[i].set_title(f'{titles[i]} Degree',fontsize=30,weight='bold',pad=15)
        axs[i].set_xlabel('Group',fontsize=25,weight='bold')
        axs[i].set_ylabel('Average Degree',fontsize=25,weight='bold')
        axs[i].spines['top'].set_linewidth(1.5)
        axs[i].spines['bottom'].set_linewidth(1.5)
        axs[i].spines['right'].set_linewidth(1.5)
        axs[i].spines['left'].set_linewidth(1.5)
        #第三个子图：平均聚类系数的柱状图
        #注意，avg_cluster应该只有一个值，如果它是一个列表且您想显示每个treat的结果，则需要调整
        if speci=='homo':
            axs[i+1].bar(range(len(avg_clusters)), avg_clusters,color=['black','red','grey','black','red','grey'])
            avgds = avg_clusters
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i+1].set_ylim(y_min, y_max)
            # max_val = max(avg_clusters) * 1.1  # 增加10%的缓冲区
            # axs[i+1].set_ylim(0, max_val)
            for j, val in enumerate(avg_clusters):
                axs[i+1].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold')

            axs[i+1].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i+1].set_xticklabels(labels,fontsize=22,weight='bold',rotation=0)
            axs[i+1].set_title(f'{titles[i]} Cluster Coef',fontsize=30,weight='bold',pad=15)
            axs[i+1].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i+1].set_ylabel('Coefficient',fontsize=25,weight='bold')
            axs[i+1].spines['top'].set_linewidth(1.5)
            axs[i+1].spines['bottom'].set_linewidth(1.5)
            axs[i+1].spines['right'].set_linewidth(1.5)
            axs[i+1].spines['left'].set_linewidth(1.5)
        else:
            axs[i+3].bar(range(len(avg_clusters)), avg_clusters,color=['black','red','grey','black','red','grey'])
            # max_val = max(avg_clusters) * 1.1  # 增加10%的缓冲区
            # axs[i+3].set_ylim(0, max_val)
            avgds = avg_clusters
            y_min = np.array(avgds).min() - (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            y_max = np.array(avgds).max() + (np.array(avgds).max() - np.array(avgds).min()) * 0.1
            axs[i+3].set_ylim(y_min, y_max)
            for j, val in enumerate(avg_clusters):
                axs[i + 3].text(j, val, f'{val:.3f}', ha='center', va='bottom',fontsize=22,weight='bold')

            axs[i+3].set_xticks(range(len(labels)))
            # 设置x轴的刻度标签
            axs[i+3].set_xticklabels(labels,fontsize=22,weight='bold',rotation=0)
            axs[i+3].set_title(f'{titles[i]} Cluster Coef',fontsize=30,weight='bold',pad=15)
            axs[i+3].set_xlabel('Group',fontsize=25,weight='bold')
            axs[i+3].set_ylabel('Coefficient',fontsize=25,weight='bold')
            axs[i+3].spines['top'].set_linewidth(1.5)
            axs[i+3].spines['bottom'].set_linewidth(1.5)
            axs[i+3].spines['right'].set_linewidth(1.5)
            axs[i+3].spines['left'].set_linewidth(1.5)
            # axs[0].set_ylim(23,25)
            # axs[1].set_ylim(27,28.5)
            # axs[2].set_ylim(12,13)
            # axs[3].set_ylim(0.6,0.7)
            # axs[4].set_ylim(0.5,0.6)
            # axs[5].set_ylim(0.3,0.4)
        i=i+1

    plt.tight_layout()
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/submodule_degree/{data_type}_allGRN_dgree3.png',bbox_inches='tight')








