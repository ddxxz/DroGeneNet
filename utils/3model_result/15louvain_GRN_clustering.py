import networkx as nx
from community import community_louvain  # 需要先安装python-louvain库
import pandas as pd
import igraph as ig
import matplotlib.pyplot as plt
from infomap import Infomap
import random
from icecream import install,ic
install()
import os

#,'CE'
all_weights=[]
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']
    treats=['CD']
    for n, speci in enumerate(species):
        avgds = []
        d_hiss = []
        avg_clusters = []
        for treat in treats:
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')
            deg_up_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_up.csv')
            deg_down_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_down.csv')
            deg_gene=pd.concat([deg_up_data['FEATURE_ID'],deg_down_data['FEATURE_ID']],axis=0)
            data = data[(data['TF'].isin(deg_gene) | data['Target'].isin(deg_gene))]

            G = nx.DiGraph()
            G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
            # if speci=='indica':
            #     resolution=1.038
            # elif speci=='japonica':
            #     resolution=1.2
            # elif speci=='zeamays':
            #     resolution=1.5
            if speci=='indica_BGI':
                resolution= 1.2   #0.6
            elif speci=='japonica':
                resolution=1.2
            elif speci=='zeamays':
                resolution= 1.2    #0.9
            #resolution=1.5
            resolution= 1.2 
            partition = community_louvain.best_partition(G.to_undirected(), resolution=resolution, randomize=False)
            ic(len(partition.keys()))

            # 收集每个社区的节点
            community_dict = {}
            for node, community_id in partition.items():
                community_dict.setdefault(community_id, []).append(node)

            #sampled_G = sample_graph(G, fraction=0.01)  # 采样 1% 的节点
            #draw_graph(G, community_dict, f"{speci} - {treat}",speci,treat)
            ic(f"\n总社区个数: {len(set(partition.values()))}")

            if treat=='CD':
                # 保存每个社区到一个 CSV 文件
                for community_id, nodes in community_dict.items():
                    subgraph_edges = []
                    for node in nodes:
                        for target in G.successors(node):
                            if target in nodes:
                                edge_weight = G[node][target]['pred_weights2']
                                subgraph_edges.append([node, target, edge_weight])
                    
                    #仅保存包含边的社区
                    if subgraph_edges:
                        df = pd.DataFrame(subgraph_edges, columns=['Source', 'Target', 'Weight'])
                        newdir=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/louvein_GRNcluster/{data_type}_{speci}/CK-CD/'
                        if not os.path.exists(newdir):
                            os.makedirs(newdir)
                        csv_filename = f"{newdir}/{speci}_{treat}_{community_id}.csv"
                        df.to_csv(csv_filename, index=False)
                        print(f"Saved community {community_id} to {csv_filename}")
                    else:
                        print(f"Community {community_id} has no edges and will not be saved.")

    #species=['indica_BGI','japonica','zeamays']
    treats=['CE']#,'CE'
    all_weights=[]

    for n, speci in enumerate(species):
        avgds = []
        d_hiss = []
        avg_clusters = []
        for treat in treats:
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')
            deg_up_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG_up.csv')
            deg_down_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG_down.csv')
            deg_gene=pd.concat([deg_up_data['FEATURE_ID'],deg_down_data['FEATURE_ID']],axis=0)
            data = data[(data['TF'].isin(deg_gene) | data['Target'].isin(deg_gene))]

            G = nx.DiGraph()
            G = nx.from_pandas_edgelist(data, 'TF', 'Target', edge_attr='pred_weights2', create_using=nx.DiGraph())
            # if speci=='indica':
            #     resolution=1.038
            # elif speci=='japonica':
            #     resolution=1.2
            # elif speci=='zeamays':
            #     resolution=1.5
            if speci=='indica_BGI':
                resolution= 1.2   #0.6
            elif speci=='japonica':
                resolution=1.2
            elif speci=='zeamays':
                resolution= 1.2    #0.9
            #resolution=1.5
            resolution= 1.2
            partition = community_louvain.best_partition(G.to_undirected(), resolution=resolution, randomize=False)
            ic(len(partition.keys()))

            # 收集每个社区的节点
            community_dict = {}
            for node, community_id in partition.items():
                community_dict.setdefault(community_id, []).append(node)

            #sampled_G = sample_graph(G, fraction=0.01)  # 采样 1% 的节点
            #draw_graph(G, community_dict, f"{speci} - {treat}",speci,treat)
            ic(f"\n总社区个数: {len(set(partition.values()))}")

            if treat=='CE':
                # 保存每个社区到一个 CSV 文件
                for community_id, nodes in community_dict.items():
                    subgraph_edges = []
                    for node in nodes:
                        for target in G.successors(node):
                            if target in nodes:
                                edge_weight = G[node][target]['pred_weights2']
                                subgraph_edges.append([node, target, edge_weight])
                    
                    #仅保存包含边的社区
                    if subgraph_edges:
                        df = pd.DataFrame(subgraph_edges, columns=['Source', 'Target', 'Weight'])
                        newdir=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/louvein_GRNcluster/{data_type}_{speci}/CK-CE/'
                        if not os.path.exists(newdir):
                            os.makedirs(newdir)
                        csv_filename = f"{newdir}/{speci}_{treat}_{community_id}.csv"
                        df.to_csv(csv_filename, index=False)
                        print(f"Saved community {community_id} to {csv_filename}")
                    else:
                        print(f"Community {community_id} has no edges and will not be saved.")
