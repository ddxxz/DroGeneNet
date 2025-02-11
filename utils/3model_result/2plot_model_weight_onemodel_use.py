import sys
from pathlib import Path
import networkx as nx

#软连接      ln -s ../../../../dxzDL ./dxzDL
sys.path.append('/home/win/4T/GeneDL/DXZ_DL/expriments/GRN')
#import shap
import numpy as np
import torch
import os
import csv
import yaml
import networkx as nx
import pandas as pd
from torch.utils.data import DataLoader,Dataset

import matplotlib.pyplot as plt
from dxzDL.models.GRN.Model_graph import Graphormer1Link,Graphormer2Link_gate,Graphormer2Link_gate_type,Graphormer1Link_type
#from models.Model_indexfind import *
from omegaconf import OmegaConf
from scipy import stats
from types import SimpleNamespace
from data_graph import *
#from main_graph_onemodel import evaluation
from icecream import install,ic
install()
from sklearn.decomposition import PCA
from sklearn.manifold import Isomap
import pickle
import umap
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
from sklearn.metrics import calinski_harabasz_score
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage
config = {
    "font.family":'Times New Roman',
    "font.size": 18,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
}
plt.rcParams.update(config)
# base_path = '/mnt/e/deep_learning/outputs/outputs20231020/att235_resample_compress_mask_gumbel_attloss/'
# model_results_path = '/mnt/e/deep_learning/outputs/outputs20231020/att235_resample_compress_mask_gumbel_attloss/model_results/Indexfind_onehot_adacos_ gumbel/base'
def embedding_plot(speci,condition,data,name):
    # 降维 - PCA
    # isomap = Isomap(n_components=2, n_neighbors=20)  # n_neighbors 控制每个点邻近点的数量
    # data_isomap = isomap.fit_transform(data)

    # 降维 - t-SNE
    tsne = TSNE(n_components=2, perplexity=50, early_exaggeration=30, random_state=42)
    data_tsne = tsne.fit_transform(data)

    # 降维 - UMAP
    reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, n_components=2, random_state=42)
    data_umap = reducer.fit_transform(data)

    # 聚类 - 使用层次聚类
    def perform_clustering(data, threshold):
        # 使用 AgglomerativeClustering 进行聚类
        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=threshold, linkage='ward')
        labels = clustering.fit_predict(data)
        return labels

    # 对不同的降维数据进行聚类

    threshold_tsne = 0.3  # 假设 t-SNE 可能需要更紧密的阈值
    threshold_umap = 0.5
    #threshold_isomap = 0.5  # 假设 Isomap 可能需要和其他方法一样的阈值
    #labels_isomap = perform_clustering(data_isomap, threshold_isomap)
    labels_tsne = perform_clustering(data_tsne, threshold_tsne)
    labels_umap = perform_clustering(data_umap, threshold_umap)

    #ch_score_isomap = calinski_harabasz_score(data_isomap, labels_isomap)
    ch_score_tsne = calinski_harabasz_score(data_tsne, labels_tsne)
    ch_score_umap = calinski_harabasz_score(data_umap, labels_umap)
    # 可视化 - PCA
    plt.figure(figsize=(12, 5))
    #plt.subplot(131)
    #plt.scatter(data_isomap[:, 0], data_isomap[:, 1], c=labels_isomap, cmap='viridis')
    # plt.title('Isomap with Clusters')
    # plt.colorbar()
    # plt.text(0.1, 0.95, f'{ch_score_isomap:.1f}', transform=plt.gca().transAxes, fontsize=20,weight='bold', verticalalignment='top')


    # 可视化 - t-SNE
    plt.subplot(121)
    plt.scatter(data_tsne[:, 0], data_tsne[:, 1], c=labels_tsne, cmap='viridis')
    plt.title('t-SNE with Clusters')
    plt.colorbar()
    plt.text(0.1, 0.95, f'{ch_score_tsne:.1f}', transform=plt.gca().transAxes, fontsize=20, weight='bold', verticalalignment='top')


    # 可视化 - UMAP
    plt.subplot(122)
    plt.scatter(data_umap[:, 0], data_umap[:, 1], c=labels_umap, cmap='viridis')
    plt.title('UMAP with Clusters')
    plt.colorbar()
    plt.text(0.1, 0.95, f'{ch_score_umap:.1f}', transform=plt.gca().transAxes, fontsize=20, weight='bold', verticalalignment='top')

    plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/graphormer2Link_alltype/model_weight_embedding/{speci}_{condition}_{name}_embedding.png',bbox_inches='tight')


    #ic(f'{speci}_{condition}_ch_score',ch_score)

def evaluation(args, model, loader,data_feature,adj,genenum,att_bias):
    device = model.device
    model.eval()
    outputs = []
    labels = []
    for data in loader:
        with torch.no_grad():
            input = data['data'].to(device=device, dtype=torch.float)
            input = input.to(torch.long)
            #genenum,
            input,adj,x0,edge_index,output,tf_embed,target_embed,x1,alpha1,x2,alpha2,x3,alpha3 = model(data_feature,adj,input,att_bias)
            #ic(output)
            if args.flag:
                output = torch.softmax(output, dim=1)
            else:
                output = torch.sigmoid(output)
            outputs.append(output.cpu().numpy())
            label = data['Label'].to(device=device, dtype=torch.float)
            labels.append(label.cpu().numpy())
    #ic(outputs.shape,labels.shape)
    outputs = np.concatenate(outputs)
    labels = np.concatenate(labels)
    y_pred = outputs
    y_true = labels
    if args.flag:
        # y_p = torch.argmax(y_pred,dim=1)
        # y_p = y_pred[:,-1]
        # y_p = y_p
        y_p = y_pred.flatten()
    else:
        y_p = y_pred
        y_p = y_p.flatten()

    y_pred_binary = [1 if p >= 0.5 else 0 for p in y_p]
    y_t = y_true.flatten().astype(int)

    AUC = roc_auc_score(y_true=y_t, y_score=y_p)
    AUPR = average_precision_score(y_true=y_t,y_score=y_p)
    AUPR_norm = AUPR/np.mean(y_t)

    ACC = metrics.accuracy_score(y_true=y_t, y_pred=y_pred_binary)

    TP = calculate_TP(y_t, y_p)
    TN = calculate_TN(y_t, y_p)
    FP = calculate_FP(y_t, y_p)
    FN = calculate_FN(y_t, y_p)
    precision = calculate_precision(y_t, y_p)
    recall = calculate_recall(y_t, y_p)
    F1_score = calculate_F1(y_t, y_pred_binary)
    # log_loss = calculate_log_loss_weighted(y_t, y_p)
    # log_loss = np.mean(log_loss)
    brier_score = calculate_brier_score(y_t, y_p)
    #mrr,hits_at_1,hits_at_3,hits_at_10 = calculate_MRR_HITS(y_t,y_p)
    #score_dict = calculate_MR_MRR_HITS(y_t,y_p)
    hits1 = calculate_hits(y_p,1)
    hits3 = calculate_hits(y_p,3)
    hits10 = calculate_hits(y_p,10)
    hits30 = calculate_hits(y_p,30)
    mr,mrr = calculate_mr_mrr(y_p)
    score_dict = {
        'AUC': round(AUC, 3),
        'AUPR': round(AUPR, 3),
        'ACC': round(ACC, 3),
        'F1_score': round(F1_score, 3),
        'brier_score': round(brier_score, 3),
        'hits@1': round(hits1, 3),
        'hits@3': round(hits3, 3),
        'hits@10': round(hits10, 3),
        'hits@30': round(hits30, 3),                        
        'MR': round(mr, 3),
        'MRR': round(mrr, 3),
    }
    return outputs,labels,score_dict


conditions=['CK','CD','CE']#
species=['indica_BGI','japonica','zeamays']#,
data_types=['TF','KEGG','PPI',]#
preds=[]
tf_embeds=[]
target_embeds=[]
x1s=[]
alpha1s=[]
x2s=[]
alpha2s=[]
x3s=[]
alpha3s=[]
for data_type in data_types:
    for speci in species:
        for i,condition in enumerate(conditions):
            model_path = f'/home/win/4T/GeneDL/data_output/outputs/batchcorrrection_GRN/exp_Graphormer2Link_gate_type_all_type_{data_type}_{speci}_{speci}'
            conf_path = f'{model_path}/.hydra/config.yaml'

            args = OmegaConf.load(conf_path)
            args = SimpleNamespace(**args)
            datas=[args.expression.CK_exp,args.expression.CD_exp,args.expression.CE_exp]#
            #print(file)
            #print("resample:",args.resample,"compres:",args.compress_num)
            plt.rcParams['font.sans-serif']=['Times New Roman']#'SimHei',
            plt.rcParams['axes.unicode_minus'] = False
            #args.resample =60
            data_input = pd.read_csv(datas[i],index_col=0)
            # if speci=='indica':
            #     repeat_cols = data_input.iloc[:, :68]
            #     data_input = pd.concat([data_input, repeat_cols], axis=1)
            # elif speci=='japonica':
            #     data_input = data_input.iloc[:,:-32]
            # else:
            #     data_input = data_input.iloc[:,:-35]

            loader = load_data(data_input)
            feature = loader.exp_data()
            #ic(feature.shape,feature[0,2])
            tf = pd.read_csv(args.expression.tf_file,index_col=0)['index'].values.astype(np.int64)
            #ic(tf.shape,tf[0])
            target = pd.read_csv(args.expression.target_file,index_col=0)['index'].values.astype(np.int64)
            #ic(target.shape,target[0])
            feature = torch.from_numpy(feature)
            tf = torch.from_numpy(tf)

            device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
            data_feature = feature.to(device)
            tf = tf.to(device)



            train_data = BasicDataset(args,args.dset.tr_csv,feature.shape[0],flag=args.flag)
            train_loader = DataLoader(train_data, batch_size=args.batch_size, shuffle=False, pin_memory=True, drop_last=False)
            train_adj = train_data.Adj_Generate(tf,flag=args.flag,loop=args.loop,adj_type=args.adj_type)
            #adj=train_adj
            train_adj = adj2saprse_tensor(train_adj)
            #
            val_data = BasicDataset(args,args.dset.cv_csv,feature.shape[0],flag=args.flag)
            val_loader = DataLoader(val_data, batch_size=args.batch_size, shuffle=False, pin_memory=True, drop_last=False)

            test_data = BasicDataset(args,args.dset.tt_csv,feature.shape[0],flag=args.flag)
            test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, pin_memory=True, drop_last=False)

            G = nx.from_numpy_array(train_adj.to_dense().numpy(), create_using=nx.DiGraph)
            V = train_adj.shape[0]
            space_encoding = np.full((V, V), -1, dtype=int)

            # 计算最短路径并填充空间编码矩阵
            length = dict(nx.all_pairs_shortest_path_length(G))
            for i, targets in length.items():
                for j, path_len in targets.items():
                    space_encoding[i, j] = path_len

            # 将 space_encoding 转换为 tensor
            space_encoding_tensor = torch.tensor(space_encoding, dtype=torch.float32)
            ic(space_encoding_tensor.max())
            # 初始化 att_bias，将 -1e9 赋值给不存在路径的位置
            att_bias = torch.zeros_like(space_encoding_tensor)
            att_bias[space_encoding_tensor == -1] = -1e9

            # 如果需要将 tensor 转移到 GPU
            device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

            #space_encoding_tensor = space_encoding_tensor.to_sparse()
            space_encoding_tensor = space_encoding_tensor.to(device)
            
            att_bias = att_bias.to(device)
            train_adj = train_adj.to(device)

            inputs=[]
            for data in train_loader:
                input = data['data'].to(device=device, dtype=torch.float32)
                inputs.append(input)
            x_train = torch.concat(inputs,dim=0)

            inputs=[]
            for data in val_loader:
                input = data['data'].to(device=device, dtype=torch.float32)
                inputs.append(input)
            x_val = torch.concat(inputs,dim=0)

            inputs=[]
            for data in test_loader:
                input = data['data'].to(device=device, dtype=torch.float32)
                inputs.append(input)
            x_test = torch.concat(inputs,dim=0)

            x_data = torch.concat([x_train,x_val,x_test],dim=0)
            x_data = x_data.to(device)
            model_name='Graphormer2Link_gate_type'
            with torch.no_grad():
                adj_dense = train_adj.to_dense()
                out_degree = torch.sum(adj_dense, dim=1).long()
                max_out_degree = torch.max(out_degree)
                in_degree = torch.sum(adj_dense, dim=0).long()
                max_in_degree = torch.max(in_degree) 
            # print(models)
            #path = base_path + '/' +'band_results/Indexfind_onehot_adacos_ gumbel/sorted_softmax/resampleplot_right' + '/' +model_name
            # 创建文件夹
            device = torch.device('cuda:0')
            #os.makedirs(path, exist_ok=True)
            model = eval(model_name)(input_dim=feature.size()[1],#
                    hidden1_dim=128,
                    hidden2_dim=64,
                    hidden3_dim=32,
                    output_dim=16,
                    flag=args.flag,
                    device=device,
                    type='dot',
                    hidden_dim=data_input.shape[1],
                    visualize=True,
                    adj=train_adj,
                    max_out_degree =max_out_degree,
                    max_in_degree =max_in_degree
                            )

            model = model.to(device)

            def load_checkpoint(model, checkpoint):
                keys_list = list(checkpoint['model_state'].keys())
                for key in keys_list:
                    if 'orig_mod.' in key:
                        deal_key = key.replace('_orig_mod.', '')
                        checkpoint['model_state'][deal_key] = checkpoint['model_state'][key]
                        del checkpoint['model_state'][key]
                model.load_state_dict(checkpoint['model_state'])

            # model = torch.compile(model,mode='max-autotune')
            checkpoint = torch.load(
                f'{model_path}/best_model_both.pt',
                map_location='cpu')
            ic(checkpoint.keys())
            # model.load_state_dict(checkpoint['model_state'])
            load_checkpoint(model,checkpoint)
            model.eval()

            maskeds=[]
            ic(x_data.shape,data_feature.shape)
            #for i in range(x_test.shape[0]):
            with torch.no_grad():
                x_data = x_data.to(device)
                x_data = x_data.to(torch.long)
                #ic(train_adj.dtype,x_test.dtype)
                input,adj,x0,edge_index,pred,tf_embed,target_embed,x1,alpha1,x2,alpha2,x3,alpha3=model(data_feature,train_adj,x_data,att_bias)
                train_outputs,train_labels,train_score_dict = evaluation(args, model, train_loader,data_feature,train_adj,None,att_bias)
                directory = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}'
                if not os.path.exists(directory):
                    os.makedirs(directory)
                np.savez(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/train_out.npz',
                        labels=(train_labels),  #
                        preds=(train_outputs),  #
                        )
                ic(train_score_dict)
                    #f"train_AUC: {train_score_dict['AUC']:.3f} train_AUPR: {train_score_dict['AUPR']:.3f} train_ACC: {train_score_dict['ACC']:.3f} train_F1_score: {train_score_dict['F1_score']:.3f} train_brier_score: {train_score_dict['brier_score']:.3f} train_mr: {train_score_dict['MR']:.3f} train_mrr: {train_score_dict['MRR']:.3f} train_hits_at_1: {train_score_dict['hits@1']:.3f} train_hits_at_3: {train_score_dict['hits@3']:.3f} train_hits_at_10: {train_score_dict['hits@10']:.3f}")
                val_outputs,val_labels,val_score_dict = evaluation(args, model, val_loader,data_feature,train_adj,None,att_bias)
                np.savez(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/val_out.npz',
                        labels=(val_labels),  #
                        preds=(val_outputs),  #
                        )
                ic(val_score_dict)
                    #f"val_AUC: {val_AUC:.3f} val_AUPR: {val_AUPR:.3f} val_ACC: {val_ACC:.3f} val_F1_score: {val_F1_score:.3f} val_brier_score: {val_brier_score:.3f}  val_mr: {val_score_dict['MR']:.3f} val_mrr: {val_score_dict['MRR']:.3f} val_hits_at_1: {val_score_dict['hits@1']:.3f} val_hits_at_3: {val_score_dict['hits@3']:.3f} val_hits_at_10: {val_score_dict['hits@10']:.3f}")

                test_outputs,test_labels,test_score_dict = evaluation(args, model, test_loader,data_feature,train_adj,None,att_bias)
                np.savez(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edge/{data_type}_{speci}_{condition}/test_out.npz',
                        labels=(test_labels),  #
                        preds=(test_outputs),  #
                        )
                # ic(test_score_dict)

                # pred =pred.cpu().numpy()
                # tf_embed = tf_embed.cpu().numpy()
                # target_embed = target_embed.cpu().numpy()
                # alpha1 =alpha1.cpu().numpy()
                # alpha2 =alpha2.cpu().numpy()
                # alpha3 =alpha3.cpu().numpy()
                # input =input.cpu().numpy()
                # x0 =x0.cpu().numpy()
                # x1 =x1.cpu().numpy()
                # x2 =x2.cpu().numpy()
                # x3 =x3.cpu().numpy()
                # adj =adj.to_dense().cpu().numpy()
                # edge_index =edge_index.cpu().numpy()
                # #ic(alpha3)
                # #ic(tf_embed.shape,target_embed.shape,x1.shape,alpha1.shape,x2.shape,alpha2.shape,x3.shape,alpha3.shape)

                # num_nodes = edge_index.max().item() + 1
                # values = np.mean(alpha1,axis=1).squeeze()
                # ic(values.shape)
                # adj1 = torch.sparse_coo_tensor(edge_index, values, (num_nodes, num_nodes))
                # values = np.mean(alpha2,axis=1).squeeze()
                # adj2 = torch.sparse_coo_tensor(edge_index, values, (num_nodes, num_nodes))
                # values = np.mean(alpha3,axis=1).squeeze()
                # adj3 = torch.sparse_coo_tensor(edge_index, values, (num_nodes, num_nodes))
                # ic(adj3.to_dense().shape)
                # print(f'{speci}_{condition}',np.mean(alpha3),np.max(alpha3),np.min(alpha3))


                # pd.DataFrame(pred).to_csv(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_pred.txt',sep='\t')
                # pd.DataFrame(x0).to_csv(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_x0.txt',sep='\t')
                # pd.DataFrame(x1).to_csv(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_x1.txt',sep='\t')
                # pd.DataFrame(x2).to_csv(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_x2.txt',sep='\t')
                # pd.DataFrame(x3).to_csv(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_x3.txt',sep='\t')
                # pd.DataFrame(adj1.to_dense()).to_pickle(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_adj1.pkl')
                # pd.DataFrame(adj2.to_dense()).to_pickle(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_adj2.pkl')
                # pd.DataFrame(adj3.to_dense()).to_pickle(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_adj3.pkl')
                # pd.DataFrame(adj).to_pickle(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/{speci}_{condition}_geneembedding_adj.pkl')
                
                # dense_adj = adj3.to_dense()
                # #column_sums = dense_adj.sum(dim=1)
                # #ic(column_sums.shape)
                # # 找出列和最大的前三个列索引
                # mask = dense_adj > 0.1
                # # 使用这个布尔矩阵过滤出大于0.1的元素，其余位置为0
                # filtered_values = torch.where(mask, dense_adj, torch.tensor(0.0))
                # # 按行求和
                # column_sums = filtered_values.sum(dim=1)

                # top_three_columns = torch.topk(column_sums, 3).indices
                # print("Indices of the top three columns by sum:", top_three_columns)

                # # 对每个顶部索引找出对应列中最大的三个值的索引
                # for index in top_three_columns:
                #     column_values = dense_adj[index,: ]
                #     top_three_values_in_column = torch.topk(column_values, 5).indices
                #     print(f"Top three indices in column {index}: {top_three_values_in_column}")

                # matrix_np = adj1.to_dense().numpy()
                # height, width = matrix_np.shape
                # block_size = 100  # 每个小矩阵的尺寸
                # # 计算可以分割出的矩阵块数量
                # num_blocks_height = height // block_size
                # num_blocks_width = width // block_size

                # # 使用双层循环遍历每一个100x100的子矩阵
                # for i in range(num_blocks_height):
                #     for j in range(num_blocks_width):
                #         # 计算当前子矩阵的起始行列索引
                #         start_row = i * block_size
                #         end_row = (i + 1) * block_size
                #         start_col = j * block_size
                #         end_col = (j + 1) * block_size
                        
                #         # 提取子矩阵
                #         scaled_matrix = matrix_np[start_row:end_row, start_col:end_col]
                #         if np.all(scaled_matrix == 0):
                #             continue
                #         # scaled_matrix = np.round(scaled_matrix).astype(np.uint8)
                #         plt.figure(figsize=(10, 10))  # 可以调整图像大小
                #         cax = plt.imshow(scaled_matrix, interpolation='none', cmap='hot')
                #         cbar = plt.colorbar(cax)
                #         #plt.imshow(matrix_np, interpolation='none', cmap='viridis')
                #         #plt.colorbar()  # 显示颜色条
                #         plt.title(f'{speci}-{condition}-adj1',fontsize=30,weight='bold',pad=10)
                #         plt.xlabel('Target',fontsize=25,weight='bold')
                #         plt.ylabel('TF',fontsize=25,weight='bold')
                #         plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/attention_weight/{speci}_{condition}_adj1_{i}_{j}.png',bbox_inches='tight')
                
                # matrix_np = adj2.to_dense().numpy()
                # height, width = matrix_np.shape
                # block_size = 100  # 每个小矩阵的尺寸
                # # 计算可以分割出的矩阵块数量
                # num_blocks_height = height // block_size
                # num_blocks_width = width // block_size

                # # 使用双层循环遍历每一个100x100的子矩阵
                # for i in range(num_blocks_height):
                #     for j in range(num_blocks_width):
                #         # 计算当前子矩阵的起始行列索引
                #         start_row = i * block_size
                #         end_row = (i + 1) * block_size
                #         start_col = j * block_size
                #         end_col = (j + 1) * block_size
                        
                #         scaled_matrix = matrix_np[start_row:end_row, start_col:end_col]
                #         if np.all(scaled_matrix == 0):
                #             continue
                #         plt.figure(figsize=(10, 10))  # 可以调整图像大小
                #         cax = plt.imshow(scaled_matrix, interpolation='none', cmap='hot')
                #         cbar = plt.colorbar(cax)
                #         #plt.colorbar()  # 显示颜色条
                #         plt.title(f'{speci}-{condition}-adj2',fontsize=30,weight='bold',pad=10)
                #         plt.xlabel('Target',fontsize=25,weight='bold')
                #         plt.ylabel('TF',fontsize=25,weight='bold')
                #         plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/model_weight/attention_weight/{speci}_{condition}_adj2_{i}_{j}.png',bbox_inches='tight')
                
                # matrix_np = adj3.to_dense().numpy()
                # height, width = matrix_np.shape
                # block_size = 100  # 每个小矩阵的尺寸
                # # 计算可以分割出的矩阵块数量
                # num_blocks_height = height // block_size
                # num_blocks_width = width // block_size

                # # 使用双层循环遍历每一个100x100的子矩阵
                # for i in range(num_blocks_height):
                #     for j in range(num_blocks_width):
                #         if speci=='japonica' and i==3 and j==48:
                #         #'zeamays' and i==2 and j==149:# indica i==3 and j==2
                            
                #         # 计算当前子矩阵的起始行列索引
                #             start_row = i * block_size
                #             end_row = (i + 1) * block_size
                #             start_col = j * block_size
                #             end_col = (j + 1) * block_size

                #             scaled_matrix = matrix_np[start_row:end_row, start_col:end_col]
                #             if np.all(scaled_matrix == 0):
                #                 continue
                #             plt.figure(figsize=(10, 10))  # 可以调整图像大小
                #             cax = plt.imshow(scaled_matrix, interpolation='none', cmap='hot')
                #             cbar = plt.colorbar(cax)
                #             #plt.colorbar()  # 显示颜色条
                #             plt.title(f'{speci}-{condition}-adj3',fontsize=30,weight='bold',pad=10)
                #             plt.xlabel('Target',fontsize=25,weight='bold')
                #             plt.ylabel('TF',fontsize=25,weight='bold')
                #             plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/graphormer2Link_alltype/model_weight_attention/{speci}_{condition}_adj3_{i}_{j}.png',bbox_inches='tight')
                #         elif speci=='zeamays' and i==2 and j==149:
                #             start_row = i * block_size
                #             end_row = (i + 1) * block_size
                #             start_col = j * block_size
                #             end_col = (j + 1) * block_size

                #             scaled_matrix = matrix_np[start_row:end_row, start_col:end_col]
                #             if np.all(scaled_matrix == 0):
                #                 continue
                #             plt.figure(figsize=(10, 10))  # 可以调整图像大小
                #             cax = plt.imshow(scaled_matrix, interpolation='none', cmap='hot')
                #             cbar = plt.colorbar(cax)
                #             #plt.colorbar()  # 显示颜色条
                #             plt.title(f'{speci}-{condition}-adj3',fontsize=30,weight='bold',pad=10)
                #             plt.xlabel('Target',fontsize=25,weight='bold')
                #             plt.ylabel('TF',fontsize=25,weight='bold')
                #             plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/graphormer2Link_alltype/model_weight_attention/{speci}_{condition}_adj3_{i}_{j}.png',bbox_inches='tight')
                #         elif speci=='indica' and i==3 and j==2:
                #             start_row = i * block_size
                #             end_row = (i + 1) * block_size
                #             start_col = j * block_size
                #             end_col = (j + 1) * block_size

                #             scaled_matrix = matrix_np[start_row:end_row, start_col:end_col]
                #             if np.all(scaled_matrix == 0):
                #                 continue
                #             plt.figure(figsize=(10, 10))  # 可以调整图像大小
                #             cax = plt.imshow(scaled_matrix, interpolation='none', cmap='hot')
                #             cbar = plt.colorbar(cax)
                #             #plt.colorbar()  # 显示颜色条
                #             plt.title(f'{speci}-{condition}-adj3',fontsize=30,weight='bold',pad=10)
                #             plt.xlabel('Target',fontsize=25,weight='bold')
                #             plt.ylabel('TF',fontsize=25,weight='bold')
                #             plt.savefig(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/graphormer2Link_alltype/model_weight_attention/{speci}_{condition}_adj3_{i}_{j}.png',bbox_inches='tight')

                
                
                
                #data_feature = data_feature.cpu().numpy()
                # embedding_plot(speci,condition,input,'input')
                # embedding_plot(speci,condition,x0,'x0')
                # embedding_plot(speci,condition,x1,'x1')
                # embedding_plot(speci,condition,x2,'x2')
                # embedding_plot(speci,condition,x3,'x3')

           

        # tf_embeds.append(tf_embed)
        # target_embeds.append(target_embed)
        # alpha1s.append(alpha1)
        # alpha2s.append(alpha2)
        # alpha3s.append(alpha3)
        # x1s.append(x1)
        # x2s.append(x2)
        # x3s.append(x3)
        #pd.DataFrame(alpha3s).to_csv(f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_gene_embedding_WGCNA/{speci}_{condition}_alpha3s.csv',index=False)

    #i.shape: torch.Size([2, 258637])
    #v.shape: torch.Size([258637])
    # tf_embed.shape: (342523, 16)
    # target_embed.shape: (342523, 16)
    # x1.shape: (27726, 384)
    # alpha1.shape: (258637, 3)
    # x2.shape: (27726, 192)
    # alpha2.shape: (258637, 3)
    # x3.shape: (27726, 16)
    # alpha3.shape: (258637, 1)
    #x1.shape: torch.Size([27726, 640]) 一共多少个节点，隐藏特征维度
    # alpha1.shape: (258637, 5)  多少条存在的边，注意力头数
    #x2.shape: torch.Size([27726, 192]) 一共多少个节点，隐藏特征维度
    # alpha2.shape: (258637, 3)  多少条存在的边，注意力头数
    #x3.shape: torch.Size([27726, 16]) 一共多少个节点，隐藏特征维度
    # alpha3.shape: (258637, 1)  多少条存在的边，注意力头数
    #tf_embed.shape: (342523, 16) 一共需要预测多少条边，每条边的特征维度
    #target_embed.shape: (342523, 16) 一共需要预测多少条边，每条边的特征维度
# x1 = np.concatenate(x1s, axis=1)
# x2 = np.concatenate(x2s, axis=1)
# x3 = np.concatenate(x3s, axis=1)

# # 对于注意力权重，我们纵向拼接
# alpha1 = np.concatenate(alpha1s, axis=0)
# alpha2 = np.concatenate(alpha2s, axis=0)
# alpha3 = np.concatenate(alpha3s, axis=0)

# tf_embeds = np.concatenate(tf_embed, axis=0)
# target_embeds = np.concatenate(target_embed, axis=0)

    # # 节点特征可视化
    # def visualize_node_features(x, title):
    #     tsne = TSNE(n_components=2, random_state=0)
    #     x_2d = tsne.fit_transform(x)
        
    #     plt.figure(figsize=(10, 6))
    #     plt.scatter(x_2d[:, 0], x_2d[:, 1], alpha=0.5)
    #     plt.title(title)
    #     plt.savefig(f'/mnt/e/GeneDL/OSDrought_GCNAT_Link/plot/outfig/model_weight/{condition}_{title}.png',bbox_inches='tight')

    # # 注意力权重可视化
    # def visualize_attention_weights(alpha, title):
    #     plt.figure(figsize=(12, 8))
    #     sns.heatmap(alpha, cmap='viridis')
    #     plt.title(title)
    #     plt.savefig(f'/mnt/e/GeneDL/OSDrought_GCNAT_Link/plot/outfig/model_weight/{condition}_{title}.png',bbox_inches='tight')

    # # 节点特征
    # visualize_node_features(x1, f"{condition}-Node-Features-Visualization-for-x1")
    # visualize_node_features(x2, f"{condition}-Node-Features-Visualization-for-x2")
    # visualize_node_features(x3, f"{condition}-Node-Features-Visualization-for-x3")

    # # 注意力权重
    # visualize_attention_weights(alpha1, f"{condition}-Attention Weights Visualization for alpha1")
    # visualize_attention_weights(alpha2, f"{condition}-Attention Weights Visualization for alpha2")
    # visualize_attention_weights(alpha3, f"{condition}-Attention Weights Visualization for alpha3")

    # # 边特征散点图比较
    # visualize_node_features(tf_embed,f"TF Embed")
    # visualize_node_features(target_embed,f"Target Embed")
    # plt.figure(figsize=(10, 6))
    # plt.scatter(tf_embed[:, 0], tf_embed[:, 1], alpha=0.5, label='TF Embed')
    # plt.scatter(target_embed[:, 0], target_embed[:, 1], alpha=0.5, label='Target Embed')
    # plt.title("Edge-Feature-Embeddings-Comparison")
    # plt.legend()
    # plt.savefig(f'/mnt/e/GeneDL/OSDrought_GCNAT_Link/plot/outfig/model_weight/{condition}_Edge-Feature-Embeddings-Comparison.png',bbox_inches='tight')






