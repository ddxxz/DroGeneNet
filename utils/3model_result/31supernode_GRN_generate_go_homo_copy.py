import pandas as pd
import ast
from collections import Counter
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.lines import Line2D
from icecream import ic
plt.figure(figsize=(10, 10))
data_types = ['PPI', 'KEGG']#'TF', 
for data_type in data_types:
    species=['homo']
    for speci in species:
        data_CD=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/CD_all.txt',sep='\t')

        deg_up_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_up.csv')
        deg_down_data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG_down.csv')
        deg_gene=pd.concat([deg_up_data['FEATURE_ID'],deg_down_data['FEATURE_ID']],axis=0)

        filtered_df = data_CD[(data_CD['TF'].isin(deg_gene) & data_CD['Target'].isin(deg_gene))]
        ic(len(filtered_df))
        #ic(filtered_df.keys())
        #ic(filtered_df['target_family'])
        category_file = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/supernode_GRN/{data_type}_{speci}_CK_CD_go_class.csv'
        category_numbers = pd.read_csv(category_file, header=None).squeeze().tolist()#[1:]
        #ic(category_numbers)
        #数据预处理
        #ic(filtered_df['tf_family']) 
        ic(len(filtered_df['TF']),len(filtered_df['TF'].unique()))
        #filtered_df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/supernode_GRN/test_{speci}.csv',index=False)
        tf_nodes = filtered_df.groupby('tf_family')['TF'].apply(list).to_dict()
        #ic(tf_nodes)
        tf_nodes = {k: v for k, v in tf_nodes.items() if k != 'target_gene'}
        top_tf_nodes = dict(sorted(tf_nodes.items(), key=lambda item: len(item[1]), reverse=True)[:12])#
        #ic(top_tf_nodes)
        #for key, nodes in tf_nodes.items():
        # print(f'{key},{len(nodes)}')
        # 处理Target节点
        target_nodes = {}
        tf_interactions = {}

        for index, row in filtered_df.iterrows():
            target = row['Target']
            target_family = row['target_family']
            if target_family == 'target_gene':
                terms = ast.literal_eval(row['Target_term'])
                for go_term in terms:
                    if go_term in target_nodes:
                        target_nodes[go_term].append(target)
                    else:
                        target_nodes[go_term] = [target]
                if row['TF'] in tf_interactions:
                    tf_interactions[row['TF']].extend([go_term for go_term in terms if go_term in target_nodes])
                else:
                    tf_interactions[row['TF']] = [go_term for go_term in terms if go_term in target_nodes]
            else:
                if target_family in tf_nodes:
                    tf_nodes[target_family].append(target)
                if row['TF'] in tf_interactions:
                    tf_interactions[row['TF']].append(target)
                else:
                    tf_interactions[row['TF']] = [target]

        new_target_nodes = {str(i): [] for i in range(1, 12)}
        target_to_category = {}

        for idx, go_term in enumerate(target_nodes.keys()):
            #ic(idx)
            #ic(len(target_nodes.keys()))
            category = category_numbers[idx]
            new_target_nodes[str(category)].extend(target_nodes[go_term])
            target_to_category[go_term] = str(category)

        G = nx.Graph()
        
        # 添加源节点
        for tf_family, tfs in top_tf_nodes.items():
            G.add_node(tf_family, size=len(tfs), node_type='tf')

        # 添加目标节点
        for category, genes in new_target_nodes.items():
            G.add_node(category, size=len(genes), node_type='target')

        # 添加边
        added_edges = set()
        all_tf = []

        for k, v in top_tf_nodes.items():
            all_tf.extend(v)
        ic(len(all_tf), len(np.unique(all_tf)))

        for index, row in filtered_df.iterrows():
            tf = row['TF']
            target = row['Target']
            
            # 找到TF和target对应的超节点
            tf_super_node = None
            target_super_node = None
            
            for tf_family, tfs in top_tf_nodes.items():
                #ic(tf)
                if tf in tfs:
                    tf_super_node = tf_family
                if target in tfs:
                    target_super_node = tf_family
            
            target_categories = []
            for go_term in ast.literal_eval(row['Target_term']):
                target_category = target_to_category.get(go_term, None)
                if target_category:
                    target_categories.append(target_category)
            
            if tf_super_node:
                for target_category in target_categories:
                    edge = (tf_super_node, target_category)
                    if edge not in added_edges:
                        G.add_edge(*edge)
                        added_edges.add(edge)

            tf_categories = [tf_super_node]#[]

                # tf_category = top_tf_nodes.get(go_term, None)
                # if tf_category:
                #     tf_categories.append(go_term)

            if tf_super_node and target_super_node:
                ic(tf_super_node,target_super_node)
                #for tf_category in tf_categories:
                edge = (tf_super_node, target_super_node)
                if edge not in added_edges:
                    #ic(edge)
                    G.add_edge(*edge)
                    added_edges.add(edge)
        #ic(all_tf)
        # 确保所有节点都有 size 和 node_type 属性
        for node in G.nodes():
            if 'size' not in G.nodes[node]:
                G.nodes[node]['size'] = 1
            if 'node_type' not in G.nodes[node]:
                G.nodes[node]['node_type'] = 'unknown'

        target_types = ['Oxidative detoxification', 'Osmotic regulation', 'Plasma membrane and cell wall', 
                        'Respiratory and metabolic enzyme activity', 'protein modification', 'DNA replication and transcription',
                        'Growth and development', 'Stress response', 'Hormone regulation', 
                        'photosynthesis', 'stomata']
        colors = ['blue', 'cyan', 'magenta', 'yellow', 'orange', 'purple', 'brown', 'pink', 'gray', 'lime', 'red']


        # 可视化超图
        #pos = nx.spring_layout(G)
        # pos = nx.spring_layout(G)
        #pos = nx.circular_layout(G)
        #pos = nx.kamada_kawai_layout(G)
        #pos = nx.spectral_layout(G)
        #pos = nx.shell_layout(G)
        #pos = nx.spiral_layout(G)

        # 收集节点类型
        source_nodes = [node for node, attr in G.nodes(data=True) if attr['node_type'] == 'tf']
        target_nodes = [node for node, attr in G.nodes(data=True) if attr['node_type'] == 'target']

        # 获取每种节点类型的大小和颜色
        source_node_sizes = [G.nodes[node]['size'] * 0.5 for node in source_nodes]
        source_node_colors = ['red' for node in source_nodes]

        target_node_sizes = [G.nodes[node]['size'] * 0.5 for node in target_nodes]
        target_node_colors = colors

        pos = nx.shell_layout(G, nlist=[source_nodes, target_nodes])

        plt.figure(figsize=(10, 10))

        # 绘制源节点（三角形）
        nx.draw_networkx_nodes(G, pos, nodelist=source_nodes, node_size=source_node_sizes, node_color=source_node_colors, node_shape='^')

        # 绘制目标节点（圆形）
        for target_type, color in zip(target_types, colors):
            #ic(color)
            nx.draw_networkx_nodes(G, pos, nodelist=target_nodes, node_size=target_node_sizes, node_color=colors, node_shape='o', label=target_type)

        # 绘制边
        nx.draw_networkx_edges(G, pos, edge_color='grey', arrows=True, arrowstyle='-|>', arrowsize=20)
        # 绘制标签
        nx.draw_networkx_labels(G, pos, font_size=20, font_weight='bold')
        legend_elements = [Line2D([0], [0], marker='o', color='w', label=target_type, markersize=10, markerfacecolor=color) 
                    for target_type, color in zip(target_types, colors)]
        plt.legend(handles=legend_elements, loc='upper left', fontsize=10, labelspacing=0.1, handletextpad=0.5)

        #plt.legend(scatterpoints=1, loc='upper left', fontsize=10, markerscale=0.1, labelspacing=0.5, handletextpad=0.5)
        plt.title(f'Super GRN for {speci}',fontsize=30,fontweight='bold')
        plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/{data_type}_{speci}_CK_CD.png', bbox_inches='tight')
        plt.clf()
