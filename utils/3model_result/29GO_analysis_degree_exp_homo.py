import gseapy as gp
import os
from plot_new import barplot, dotplot
import pandas as pd
import matplotlib.pyplot as plt
from icecream import install,ic
install()
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['mathtext.fontset'] = 'stix'#和Times new roman 最像

#-------------------------------获取Geneid和GOterm---------------------------------------------------
# 假设文件是以制表符分隔的，将 'sep' 参数替换成 ',' 如果是逗号分隔的CSV
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    allspecies=['homo']
    for species in allspecies:#'homo'
    #species='zeamays'
    #species='japonica'
        if species=='zeamays':
            file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zma_GO_annotation'
        elif species=='indica_BGI':
            file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_Osi_GO_annotation'
        else:
            file_path ='/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/GOcluster/Osj_GO_annotation' # 或者 'Osj_GO_annotation.csv'

        # 读取文件，只获取前两列
        # header参数假设文件有标题行，如果没有请设定为None
        GO_term_data = pd.read_csv(file_path, sep='\t', usecols=[0, 2,3], header=0)
        if species=='zeamays':
            relation = []
            # 使用open函数和with语句读取文件
            with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt", 'r') as file:
                for line in file:
                    # 分割每行的数据，只取前两个值（原始基因名和第一个映射基因名）
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:  # 确保有足够的部分读取
                        original_gene = parts[0]
                        mapped_gene = parts[1]
                        # 保存结果
                        relation.append([original_gene, mapped_gene])

            # 将列表转换为DataFrame
            df3 = pd.DataFrame(relation, columns=['Original Gene', 'Mapped Gene'])
            # 创建映射字典
            gene_id_map = pd.Series(df3['Mapped Gene'].values, index=df3['Original Gene'].values).to_dict()
            #ic(gene_id_map)
            #GO_term_data['zeaGeneid'] = GO_term_data['Gene_id'].apply(convert_msu_to_rap)
            # 替换 column1, column2, column3 中的基因名
            GO_term_data['zeaGeneid'] = GO_term_data['Gene_id'].map(gene_id_map).fillna(GO_term_data['Gene_id'])  # fillna用于保留未找到映射的原始值
            # relation={}

            def create_gene_sets(df):
                """
                基于DataFrame创建基因集字典。DataFrame应包含两列：GO_term和OSGeneid。
                :param df: 包含GO_term和OSGeneid的DataFrame。
                :return: 基因集字典，键是GO_term，值是OSGeneid列表。
                """
                gene_sets = {}
                grouped = df.groupby('GO_term')
                
                for name, group in grouped:
                    gene_sets[name] = group['zeaGeneid'].tolist()
                return gene_sets
            
        elif species=='indica_BGI':
            GO_term_data['OSGeneid'] = GO_term_data['Gene_id']
            def create_gene_sets(df):
                """
                基于DataFrame创建基因集字典。DataFrame应包含两列：GO_term和OSGeneid。
                
                :param df: 包含GO_term和OSGeneid的DataFrame。
                :return: 基因集字典，键是GO_term，值是OSGeneid列表。
                """
                gene_sets = {}
                grouped = df.groupby('GO_term')
                for name, group in grouped:
                    gene_sets[name] = group['OSGeneid'].tolist()
                return gene_sets
            
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

            # 转换函数: 如果找不到对应的MSU标识符，则返回原RAP标识符
            def convert_rap_to_msu(rap_id):
                return relation.get(rap_id, rap_id)
            def convert_msu_to_rap(msu_id):
                return relation.get(msu_id, msu_id)

            GO_term_data['OSGeneid'] = GO_term_data['Gene_id'].apply(convert_rap_to_msu)

            def create_gene_sets(df):
                """
                基于DataFrame创建基因集字典。DataFrame应包含两列：GO_term和OSGeneid。
                
                :param df: 包含GO_term和OSGeneid的DataFrame。
                :return: 基因集字典，键是GO_term，值是OSGeneid列表。
                """
                gene_sets = {}
                grouped = df.groupby('GO_term')
                
                for name, group in grouped:
                    gene_sets[name] = group['OSGeneid'].tolist()
                
                return gene_sets
        #ic(len(GO_term_data['Gene_id'].unique()))
        ic(GO_term_data)
        #ic(GO_term_data[GO_term_data['zeaGeneid'].notna()]['zeaGeneid'])
        grouped_dataframes = {aspect: group for aspect, group in GO_term_data.groupby('Aspect')}
        # GO_term_data = GO_term_data[GO_term_data['Aspect'] == 'P']
        # all_gene_sets=create_gene_sets(GO_term_data)

        all_gene_sets=[]
        for aspect, sub_df in grouped_dataframes.items():
            #ic(f"Aspect: {aspect}")
            #ic(sub_df)

            gene_sets=create_gene_sets(sub_df)
            for key in gene_sets.keys():
                gene_sets[key] = [item for item in gene_sets[key] if item != 'None']
            all_gene_sets.append(gene_sets)

        data_path=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/'
        files=os.listdir(data_path)

        funcs=['CC','MF','BP']#
        #labels=['up','down']
        treats=['ck_cd','ck_ce']
        colors=['#2e86de','#c0392b','#28b463']
    # treats=['leaf','root']#['CK','CD']#
        #for treat in treats:
            #for label in labels:
        #     for file in files:
        all_results=[]
        #for func in funcs:
        for i,gene_sets in enumerate(all_gene_sets):
            for treat in treats:
            #for label in labels:
            #for file in files:
                #data = pd.read_csv(f'{data_path}{file}')
                #ic(data)
                #ic(gene_sets)
                # data1=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/homo/EXP_DEG/homo_CK_CD_DEG_up.csv')
                # data2=pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/new/homo/EXP_DEG/homo_CK_CD_DEG_down.csv')
                data=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/{data_type}_{species}_{treat}_degree_exp_intersection.csv')
                gene_list = data['0']
                #data = data[data['sig'] == 'sig']
                #ic(data)
                #gene_list=pd.concat([data1['FEATURE_ID'],data2['FEATURE_ID']],axis=0)
                #gene_list = pd.concat([data['Source'],data['Target']],axis=0)
                gene_list = [str(i) for i in gene_list]
                #gene_list = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/data/top_diff_CD_CE_degree.csv')
                #gene_list = list(gene_list)
                #gene_list = gene_list['0'].apply(convert_msu_to_rap)
                #ic(gene_list)
                # ic(gene_sets)
                if len(gene_list) < 5:
                    continue
                # GO富集
                #title_name=file.split('.')[0]
                enr = gp.enrichr(gene_list=gene_list,#所需查询gene_list，可以是一个列表，也可为文件（一列，每行一个基因）
                                gene_sets= gene_sets,#['GO_Biological_Process_2018'],#,#gene set library，多个相关的gene set 。如所有GO term组成一个gene set library.
                                #organism='rice',#持(human, mouse, yeast, fly, fish, worm)， 自定义gene_set 则无影响。
                                #description='kegg',#工作运行描述
                                outdir=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_GO/{data_type}_{species}',#输出目录
                                top_term=20,
                                cutoff=0.9,#pvalue阈值
                                )
                #ic(enr.results)

                #all_results.append(enr.results)
                dot_png = f"exp_{species}_{funcs[i]}_{treat}" + "_" + "dot" + ".png"
                #bar_png = f"exp_{species}_{func}" + "_" + "bar" + ".png"
                base_path = f"/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_GO/funcGO"

                if os.path.exists(base_path+dot_png):
                    os.remove(base_path+dot_png)
                #module_name=title_name.split('_')[-1]
                #fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 10))
                cut=0.8
                df= dotplot(enr.results,
                        column="Adjusted P-value",
                        size=5,
                        #figsize=(3,5),
                        xticklabels_rot=45,
                        title=f'GO Top15 Pathway DEG {species} {funcs[i]} {treat}', 
                        show_ring=True,
                        #marker='o',
                        cmap='RdBu',
                        top_term=15, legend="r",pval=0.9,
                        ofname=base_path+dot_png,
                        cutoff=cut,
                        )
                # if os.path.exists(base_path+bar_png):
                #     os.remove(base_path+bar_png)

                # ax = barplot(enr.results,
                #             column="Adjusted P-value",
                #             title=f"{species} {treat} module{module_name}",
                #             #group='Gene_set', # set group, so you could do a multi-sample/library comparsion
                #             size=1,
                #             top_term=10,
                #             #figsize=(3,5),
                #             color=colors[i],#['#27ae60', '#34495e'] , # set colors for group
                #             ofname=base_path+bar_png,
                #             cutoff=cut
                #             #color = {'KEGG_2021_Human': 'salmon', 'MSigDB_Hallmark_2020':'darkblue'}
                #             )
                df = df.iloc[::-1].reset_index(drop=True)
                df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_GO/{data_type}_{species}_{funcs[i]}_{treat}_GO.csv',index=False)

    # funcs=['CC','MF','BP']
    # label=['up','down']
    # colors=['#2e86de']#,'#c0392b','#28b463'
    # #for cla in label:
    # for data_name in data_list:
    #     all_results=[]
    #     #for i,gene_sets in enumerate(all_gene_sets):
    #     data1 = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/GRN_DEG/{species}_{data_name}_DEG_up.csv')#['DAB2IP', 'AKT1', 'ARNT', 'FLCN', 'XRCC6', 'BTG1', 'STAT1', 'STAT1']
    #     data2 = pd.read_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/GRN_DEG/{species}_{data_name}_DEG_down.csv')#['DAB2IP', 'AKT1', 'ARNT', 'FLCN', 'XRCC6', 'BTG1', 'STAT1', 'STAT1']
    #     data=pd.concat([data1,data2],axis=0)
    #     ic(data)
    #     #data = data[data['sig'] == 'sig']
    #     #ic(data)
    #     gene_list = data['Unnamed: 0']#[:5000]##['GeneID']
    #     gene_list = [str(i) for i in gene_list]
    #     #gene_list = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/data/top_diff_CD_CE_degree.csv')
    #     #gene_list = list(gene_list)
    #     #gene_list = gene_list['0'].apply(convert_msu_to_rap)
    #     # ic(gene_list)
    #     # ic(gene_sets)
    #     # GO富集
    #     enr = gp.enrichr(gene_list=gene_list,#所需查询gene_list，可以是一个列表，也可为文件（一列，每行一个基因）
    #                     gene_sets= all_gene_sets,#['GO_Biological_Process_2018'],#,#gene set library，多个相关的gene set 。如所有GO term组成一个gene set library.
    #                     #organism='rice',#持(human, mouse, yeast, fly, fish, worm)， 自定义gene_set 则无影响。
    #                     #description='kegg',#工作运行描述
    #                     outdir=f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/GO/{species}_{data_name}',#输出目录
    #                     top_term=20,
    #                     cutoff=0.9,#pvalue阈值
    #                     )
    #     ic(enr.results)
    #     #all_results.append(enr.results)
    #     dot_png = f"{species}_{data_name}" + "_" + "dot" + ".png"
    #     bar_png = f"{species}_{data_name}" + "_" + "bar" + ".png"
    #     base_path = f"/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/GO/funcGO"

    #     if os.path.exists(base_path+dot_png):
    #         os.remove(base_path+dot_png)

    #     #fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(5, 10))
    #     cut=0.8
    #     df= dotplot(enr.results,
    #             column="Adjusted P-value",
    #             size=5,
    #             #figsize=(3,5),
    #             xticklabels_rot=45,
    #             title=f'GO Top10 Pathway {species} {data_name}', 
    #             show_ring=True,
    #             #marker='o',
    #             cmap='RdBu',
    #             top_term=15, legend="r",pval=0.9,
    #             ofname=base_path+dot_png,
    #             cutoff=cut,
    #             )
    #     if os.path.exists(base_path+bar_png):
    #         os.remove(base_path+bar_png)

    #     ax = barplot(enr.results,
    #                 column="Adjusted P-value",
    #                 title=f"{species} {data_name}",
    #                 #group='Gene_set', # set group, so you could do a multi-sample/library comparsion
    #                 size=1,
    #                 top_term=20,
    #                 #figsize=(3,5),
    #                 color=colors,#['#27ae60', '#34495e'] , # set colors for group
    #                 ofname=base_path+bar_png,
    #                 cutoff=cut
    #                 #color = {'KEGG_2021_Human': 'salmon', 'MSigDB_Hallmark_2020':'darkblue'}
    #                 )
    #     df = df.iloc[::-1].reset_index(drop=True)
    #     df.to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_modelresult/{species}/GO/{species}_{data_name}.csv',index=False)


        # from gseapy import enrichment_map
        # import numpy as np
        # # return two dataframe
        # ic(enr.results)
        # nodes, edges = enrichment_map(enr.results,cutoff=cut,top_term=20)
        # ic(edges)
        # import networkx as nx
        # # build graph
        # G = nx.from_pandas_edgelist(edges,
        #                             source='src_idx',
        #                             target='targ_idx',
        #                             edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

        # fig, ax = plt.subplots()

        # # init node cooridnates
        # pos=nx.layout.spiral_layout(G)
        # ic(pos.keys())
        # #node_size = nx.get_node_attributes()
        # # draw node
        # nodelist = list(G)
        # xy = np.asarray([pos[v] for v in nodelist])
        # ic(nodelist)
        # original_list = list(xy[:, 0])
        # nx.draw_networkx_nodes(G,
        #                     pos=pos,
        #                     cmap=plt.cm.RdYlBu,
        #                     node_color= [number * 100 for number in original_list],
        #                     node_size= [number * 1000 for number in original_list])
        # # draw node label
        # ic(nodes.Term)
        # ic(pos)
        # node_term_dict = nodes.Term.to_dict()
        # labels = {v:node_term_dict[v] for v in nodelist}
        # nx.draw_networkx_labels(G,
        #                         pos=pos,
        #                         labels=labels)
        # # draw edge
        # edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
        # nx.draw_networkx_edges(G,
        #                     pos=pos,
        #                     width=list(map(lambda x: x*10, edge_weight)),
        #                     edge_color='#CDDBD4')
        # #plt.show()
        # plt.savefig(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_enrichr/GO/{species}_{data_name}_{cla}_grn.png')


# ax = dotplot(enr.res2d, title='KEGG_2021_Human',cmap='viridis_r', size=10, figsize=(3,5))
# ax = barplot(enr.res2d,title='KEGG_2021_Human', figsize=(4, 5), color='darkred')
#barplot(enr.res2d, title='CK_CD_degree',top_term=20,ofname=base_path+bar_png,)