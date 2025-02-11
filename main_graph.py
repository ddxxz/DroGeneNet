import pandas as pd
import pickle

import torch
from tqdm import tqdm
import os
import hydra
import logging

import random
import igraph as ig
from transformers import BertForMaskedLM,BertModel
from torch import optim
#from models.Model_1dref_attention import *
from models.Model_graph import *
from models.Model_demo import TRGCNGATLink,SCGCNGATLink,CAGCNGATLink,MAGLink
from models.Model_graph_foundation_downstream import GraphSAGELink_downstream,GraphSAGELink_downstream_double
from torch.utils.data import DataLoader
import torch.nn.functional as F
from torch.optim import Adam
from torch.optim.lr_scheduler import StepLR
import scipy.sparse as sp
from data_proceed.data_graph import *
from torch.utils.tensorboard import SummaryWriter
#from PytorchTools import EarlyStopping
import numpy as np
import random
import glob
import os
import time
import argparse
from icecream import install , ic
install()
from tensorboardX import SummaryWriter
from torch.utils.data import DataLoader
from sklearn.metrics import r2_score, mean_squared_error
# import torch._dynamo
# torch._dynamo.config.suppress_errors = True
import networkx as nx
logger = logging.getLogger(__name__)


def seed_torch(seed,faster=False):
    random.seed(seed)
    os.environ['PYTHONHASHSEED'] = str(seed)  # 为了禁止hash随机化，使得实验可复现
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)  # if you are using multi-GPU.

    if not faster:
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
    else:
        torch.backends.cudnn.benchmark = True
        torch.backends.cudnn.deterministic = False
seed_torch(seed=1000,faster=True)

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
            output = model(data_feature,adj,input,att_bias)
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

def embed2file(tf_embed,tg_embed,gene_file,tf_path,target_path):
    tf_embed = tf_embed.cpu().detach().numpy()
    tg_embed = tg_embed.cpu().detach().numpy()

    gene_set = pd.read_csv(gene_file, index_col=0)

    tf_embed = pd.DataFrame(tf_embed,index=gene_set['Gene'].values)
    tg_embed = pd.DataFrame(tg_embed, index=gene_set['Gene'].values)

    tf_embed.to_csv(tf_path)
    tg_embed.to_csv(target_path)

def _main(args,_run=None):
    global __file__
    # Updating paths in config
    for key, value in args.dset.items():
        if isinstance(value, str) and key not in ["matching"]:
            args.dset[key] = hydra.utils.to_absolute_path(value)
    __file__ = hydra.utils.to_absolute_path(__file__)
    tb = SummaryWriter("tensorboard")

    logger.info("For logs, checkpoints and samples check %s", os.getcwd())
    use_cuda = torch.cuda.is_available()
    device = torch.device(args.gpu_device if use_cuda else "cpu")
    print(device)

    from data_proceed.data_graph import BasicDataset
    #ic(args.expression)
    if args.data_type=='scRNA':
        ic('ok')
        data_input_dirs=args.expression.exp_file
        with open(data_input_dirs, 'rb') as f:
            data_input = pickle.load(f)
    elif args.data_type=='bulkscRNA':
        data_input_dirs=args.expression.exp_file
        with open(data_input_dirs, 'rb') as f:
            data_input_sc = pickle.load(f)
        bulk_data = pd.read_csv(args.expression.bulk_exp_file,index_col=0)
        ic(data_input_sc.shape,bulk_data.shape)
        data_input = pd.concat([data_input_sc,bulk_data],axis=1)
        #所有的
        # data_input_exp = []
        # for file in os.listdir(data_input_dirs):
        #     scdata_path = data_input_dirs + file
        #     with open(scdata_path, 'rb') as f:
        #         expression_data_df = pickle.load(f)
        #     data_input_exp.append(expression_data_df)
        # data_input = pd.concat(data_input_exp,axis=1)
        #ic(data_input.shape)
    else:
        if args.water_class == 'all':
            
            with open(args.expression.exp_file, 'rb') as file:
                # 使用 pickle 加载数据
                ic(args.expression.exp_file)
            
                data_input = pd.read_pickle(file)#pickle.load(file)
                data_input = data_input.fillna(0)
                #ic(data_input)         
            gene_name = data_input['FEATURE_ID']#['GeneID']
            data_input = data_input.iloc[:,1:]
            if args.species=='zeamays':
                gene_num=None
            else:
                with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/Geneformer/pretrain/CRR_sc_total_gene_token_dict.pickle", "rb") as fp:
                    token_dictionary = pickle.load(fp)
                gene_num = [token_dictionary[f'gene:{i}'] for i in gene_name]
                gene_num = torch.tensor(gene_num).to(device)
            ic(gene_num)
            #data_input = pd.read_csv(args.expression.exp_file,index_col=0)
        elif args.water_class == 'CD':
            data_input = pd.read_csv(args.expression.CD_exp)#,index_col=0
            gene_name = data_input['FEATURE_ID']#['GeneID']
            data_input = data_input.iloc[:,1:]
            if args.species=='zeamays':
                gene_num=None
            else:
                with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/Geneformer/pretrain/CRR_sc_total_gene_token_dict.pickle", "rb") as fp:
                    token_dictionary = pickle.load(fp)
                gene_num = [token_dictionary[f'gene:{i}'] for i in gene_name]
                gene_num = torch.tensor(gene_num).to(device)
        elif args.water_class == 'CK':
            data_input = pd.read_csv(args.expression.CK_exp) # ,index_col=0    
            gene_name = data_input['FEATURE_ID']#['GeneID']
            data_input = data_input.iloc[:,1:]
            if args.species=='zeamays':
                gene_num=None
            else:
                with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/Geneformer/pretrain/CRR_sc_total_gene_token_dict.pickle", "rb") as fp:
                    token_dictionary = pickle.load(fp)
                gene_num = [token_dictionary[f'gene:{i}'] for i in gene_name]
                gene_num = torch.tensor(gene_num).to(device)
        elif args.water_class == 'CE':
            data_input = pd.read_csv(args.expression.CE_exp)  #,index_col=0
            gene_name = data_input['FEATURE_ID']#['GeneID']
            data_input = data_input.iloc[:,1:]
            if args.species=='zeamays':
                gene_num=None
            else:
                with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/Geneformer/pretrain/CRR_sc_total_gene_token_dict.pickle", "rb") as fp:
                    token_dictionary = pickle.load(fp)
                gene_num = [token_dictionary[f'gene:{i}'] for i in gene_name]
                gene_num = torch.tensor(gene_num).to(device)


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
    train_adj = train_data.Adj_Generate(tf,args.adj_type,flag=args.flag,loop=args.loop)
    adj=train_adj
    train_adj = adj2saprse_tensor(train_adj)

    # def create_igraph_from_adj(adj_matrix):
    #     # 将稀疏矩阵转换为密集矩阵
    #     dense_adj_matrix = adj_matrix.to_dense()
    #     # 获取非零元素的索引
    #     sources, targets = dense_adj_matrix.nonzero(as_tuple=True)
    #     # 获取对应的权重
    #     weights = dense_adj_matrix[sources, targets]
    #     # 创建有向图
    #     g = ig.Graph(edges=list(zip(sources, targets)), directed=True)
    #     g.es['weight'] = weights
    #     return g

    # # 使用 igraph 计算最短路径
    # def compute_shortest_paths(g, V):
    #     space_encoding = np.full((V, V), -1, dtype=int)
    #     for i in range(V):
    #         lengths = g.shortest_paths(source=i, weights='weight', mode=ig.OUT)
    #         for j, length in enumerate(lengths[0]):
    #             if length != float('inf'):
    #                 space_encoding[i, j] = length
    #     return space_encoding
    # V = train_adj.shape[0]
    # g = create_igraph_from_adj(train_adj)
    # space_encoding = compute_shortest_paths(g, V)
    # space_encoding_tensor = torch.tensor(space_encoding, dtype=torch.float32)
    # att_bias = torch.full_like(space_encoding_tensor, 0.0)
    # att_bias[space_encoding_tensor == -1] = -1e9

    # # 如果需要将 tensor 转移到 GPU
    # device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # space_encoding_tensor = space_encoding_tensor.to(device)
    # att_bias = att_bias.to(device)
    # 使用相同的方法初始化Graph和空间编码矩阵
    # 假定 adj 是 numpy 数组格式的邻接矩阵
    G = nx.from_numpy_array(train_adj.to_dense().numpy(), create_using=nx.DiGraph)
    V = adj.shape[0]
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
    
    adj_dense = train_adj.to_dense()
    out_degree = torch.sum(adj_dense, dim=1).long()
    max_out_degree = torch.max(out_degree)
    in_degree = torch.sum(adj_dense, dim=0).long()
    max_in_degree = torch.max(in_degree) 

    # 如果需要将 tensor 转移到 GPU
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    #space_encoding_tensor = space_encoding_tensor.to_sparse()
    space_encoding_tensor = space_encoding_tensor.to(device)
    
    att_bias = att_bias.to(device)
    #ic(space_encoding_tensor.shape,att_bias.shape)
    

    val_data = BasicDataset(args,args.dset.cv_csv,feature.shape[0],flag=args.flag)
    val_loader = DataLoader(val_data, batch_size=args.batch_size, shuffle=False, pin_memory=True, drop_last=False)

    test_data = BasicDataset(args,args.dset.tt_csv,feature.shape[0],flag=args.flag)
    test_loader = DataLoader(test_data, batch_size=args.batch_size, shuffle=False, pin_memory=True, drop_last=False)


    MODELS = {
        'GCNGATLink': GCNGATLink,
        'TRGCNGATLink': TRGCNGATLink,
        'SCGCNGATLink': SCGCNGATLink,
        'CAGCNGATLink': CAGCNGATLink,
        'MAGLink': MAGLink,
        'GENELink': GENELink,
        'GNNLink': GNNLink,
        'NeighborGCNGATLink': NeighborGCNGATLink,
        'LargeGCNGATLink': LargeGCNGATLink,
        'GraphSAGELink': GraphSAGELink,
        'PinSageLink': PinSageLink,
        'Graphormer2Link_gate_type': Graphormer2Link_gate_type,
        'FastGCNLink': FastGCNLink,
        #'GraphTransformerUNet': GraphTransformerUNet,
#        'LADIESLink': LADIESLink,
        'ClusterGCNLink': ClusterGCNLink,
        'Graphormer1Link': Graphormer1Link,
        #'Graphormer1Link_gate': Graphormer1Link_gate,
        'Graphormer2Link': Graphormer2Link,
        'Graphormer1Link_gate': Graphormer1Link_gate,
        'Graphormer2Link_gate': Graphormer2Link_gate,
        'GraphSAGELink_downstream_double': GraphSAGELink_downstream_double,
        'GraphSAGELink_downstream': GraphSAGELink_downstream,
        'GraphTransformerLink': GraphTransformerLink,
        'Graphormer1Link_alphafold': Graphormer1Link_alphafold
              }
    if args.model == 'GraphSAGELink_downstream' or args.model == 'GraphSAGELink_downstream_double':                                                                  
        model = MODELS[args.model](input_dim=256,#
                    hidden1_dim=args.hidden_dim[0],
                    hidden2_dim=args.hidden_dim[1],
                    hidden3_dim=args.hidden_dim[2],
                    output_dim=args.output_dim,
                    training = args.training,
                    #num_head1=args.num_head[0],
                    num_head2=args.num_head[1],
                    flag=args.flag,
                    device=device,
                    type=args.type,
                    hidden_dim=data_input.shape[1],
                    visualize=args.visualize,
                    reduction=args.reduction,
                    adj=train_adj,
                    input_dim1=feature.size()[1],
                    )
    else:
        model = MODELS[args.model](input_dim=feature.size()[1],#
                hidden1_dim=args.hidden_dim[0],
                hidden2_dim=args.hidden_dim[1],
                hidden3_dim=args.hidden_dim[2],
                output_dim=args.output_dim,
                #training = args.training,
                #num_head1=args.num_head[0],
                #num_head2=args.num_head[1],
                flag=args.flag,
                device=device,
                type=args.type,
                hidden_dim=data_input.shape[1],
                visualize=args.visualize,
                max_out_degree=max_out_degree,
                max_in_degree=max_in_degree,
                reduction=args.reduction,
                adj=train_adj
                )
    #criterion = F.binary_cross_entropy()
    model.to(device=device)
    train_adj = train_adj.to(device)
    #train_data = train_data.to(device)
    #val_data = val_data.to(device)
    #test_data = test_data.to(device)
    logger.info(model)

    # tb.add_graph(model.cpu(), torch.zeros([1, 204]))
    if args.optim == 'Adam':
        optimizer = optim.Adam(model.parameters(),weight_decay=args.weight_decay)
    elif args.optim == 'Adadelta':
        optimizer = optim.Adadelta(model.parameters())
    elif args.optim == 'Adamax':
        optimizer = optim.Adamax(model.parameters())
    elif args.optim == 'AdamW':
        optimizer = optim.AdamW(model.parameters())
    elif args.optim == 'ASGD':
        optimizer = optim.ASGD(model.parameters())
    elif args.optim == 'LBFGS':
        optimizer = optim.LBFGS(model.parameters())
    elif args.optim == 'NAdam':
        optimizer = optim.NAdam(model.parameters(),weight_decay=args.weight_decay)
    elif args.optim == 'RAdam':
        optimizer = optim.RAdam(model.parameters())
    elif args.optim == 'RMSprop':
        optimizer = optim.RMSprop(model.parameters())
    elif args.optim == 'Rprop':
        optimizer = optim.Rprop(model.parameters())
    elif args.optim == 'SGD':
        optimizer = optim.SGD(model.parameters(),lr=0.05)
    elif args.optim == 'SparseAdam':
        optimizer = optim.SparseAdam(model.parameters())
    
    # lr = np.random.uniform(1e-3, 1e-4)
    # weight_decay = np.random.uniform(1e-2, 1e-6)
    # lr = 0
    # weight_decay = 0
    # optimizer = optim.AdamW(model.parameters())
    # optimizer = optim.RAdam(model.parameters(),weight_decay=args.weight_decay, 
    #                         lr=args.lr, betas=(0.9, 0.999), eps=1e-8)
    if args.lr_scheduler == 'StepLR':
        scheduler = optim.lr_scheduler.StepLR(optimizer,step_size=3,gamma=0.1)
    elif args.lr_scheduler == 'MultiStepLR':
        scheduler = optim.lr_scheduler.MultiStepLR(optimizer,milestones=[3,7],gamma=0.1)
    elif args.lr_scheduler == 'ExponentialLR':
        scheduler = optim.lr_scheduler.ExponentialLR(optimizer,gamma=0.1)
    elif args.lr_scheduler == 'CosineAnnealingLR':
        scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer,T_max=20)
    elif args.lr_scheduler == 'CyclicLR':
        scheduler = optim.lr_scheduler.CyclicLR(optimizer, base_lr=1e-7, max_lr=0.001,
                                            step_size_up=2000, step_size_down=2000, mode='triangular',
                                            gamma=1.0, scale_fn=None, scale_mode='cycle', cycle_momentum=False,
                                            base_momentum=0.8, max_momentum=0.9, last_epoch=-1)
    elif args.lr_scheduler == 'OneCycleLR':
        scheduler = optim.lr_scheduler.OneCycleLR(optimizer, max_lr=0.9, total_steps=100, verbose=True)
    elif args.lr_scheduler == 'CosineAnnealingWarmRestarts':
        scheduler = optim.lr_scheduler.CosineAnnealingWarmRestarts(optimizer, T_0=5, T_mult=2)
    elif args.lr_scheduler == 'ReduceLROnPlateau':
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',factor=0.1,patience=2)
    epoch = 0
    print(args.optim,args.lr_scheduler)
    if args.fine_tune:
        logger.info('-----------------------continue train------------------------')
        # torch.save(
        #     {"state": model.state_dict(),
        #      "epoch": epoch
        #      }
        #     , 'best_model.pt')
        checkpoint = torch.load(
            f'/mnt/e/deep_learning/outputs/outputs20230516/single_epoch2000weight0/contrast_Vcmax_A/{args.pretrain_name}/best_model_both.pt',
            map_location=device)
        model.load_state_dict(checkpoint['state'])
        # for p in model.parameters():#模型参数不更新
        #     p.requires_grad = False
        # for m in model.modules():#选定的bn模块参数更新,不固定
        #     if isinstance(m, nn.BatchNorm1d):
        #         m.required_grad = True

        # model.convblock1 = nn.Sequential(
        #     nn.Conv1d(in_channels=1, out_channels=50, kernel_size=5, dilation=1),
        #     nn.BatchNorm1d(num_features=50),
        #     nn.ReLU()
        # )
        # model.convblock2 = nn.Sequential(
        #     nn.Conv1d(in_channels=50, out_channels=50, kernel_size=5, dilation=2),
        #     nn.BatchNorm1d(num_features=50),
        #     nn.ReLU())
        # model.attention = nn.TransformerEncoder(
        #     nn.TransformerEncoderLayer(d_model=50, nhead=5, dim_feedforward=50 * 4, batch_first=True), num_layers=2)
        # model.fc1 = nn.Linear(in_features=400, out_features=1000)
        # # model.relu2 = nn.ReLU()
        # # model.dropout = nn.Dropout(0.2)
        # model.fc2 = nn.Linear(in_features=1000, out_features=2)
        model.to(device)

        train_r2,  train_outputs, train_labels = evalution(args, model, train_loader)
        np.savez('train_out.npz',
                 labels=(train_labels * train_data.label_std + train_data.label_mean),  #
                 preds=(train_outputs * train_data.label_std + train_data.label_mean),  #

                 )
        logger.info(f"train_direct_transfer_r2: {train_r2} ")

        val_r2, val_outputs, val_labels = evalution(args, model, val_loader)
        np.savez('val_direct_transfer_out.npz',
                 labels=(val_labels * train_data.label_std + train_data.label_mean),  #
                 preds=(val_outputs * train_data.label_std + train_data.label_mean),  #
                 )
        logger.info(f"val_direct_transfer_r2:{val_r2}")  #
        test_r2, test_outputs, test_labels = evalution(args, model, test_loader)
        np.savez('test_direct_transfer_out.npz',
                 labels=(test_labels * train_data.label_std + train_data.label_mean),  #
                 preds=(test_outputs * train_data.label_std + train_data.label_mean),  #
                 )

        logger.info(f" test_direct_transfer_r2:{test_r2}")

    elif os.path.exists('best_model_both.pt'):
        logger.info('-----------------------continue train------------------------')
        # torch.save(
        #     {"state": model.state_dict(),
        #      "epoch": epoch
        #      }
        #     , 'best_model.pt')
        checkpoint = torch.load('best_model_both.pt', map_location=device)
        model.load_state_dict(checkpoint['state'])
        epoch = checkpoint['epoch']

    def get_params(model,weight_decay = 0.):
        decay_weights, not_decay_weights = [],[]
        for name,param in model.named_parameters():
            if 'att' in name or 'ref' in name:
                not_decay_weights += [param]
            else:
                decay_weights += [param]
        params = [{"params":decay_weights,"weight_decay":weight_decay},
                  {"params": not_decay_weights, "weight_decay": 0},
                  ]
        return params
    
    # optimizer = optim.Adam(get_params(model,weight_decay=args.weight_decay), lr=args.lr, betas=(0.9, 0.999), eps=1e-8)
    # scheduler = optim.lr_scheduler.CyclicLR(optimizer, base_lr=1e-7, max_lr=0.001,
    #                                         step_size_up=2000, step_size_down=2000, mode='triangular',
    #                                         gamma=1.0, scale_fn=None, scale_mode='cycle', cycle_momentum=False,
    #                                         base_momentum=0.8, max_momentum=0.9, last_epoch=-1)

    #torch.set_float32_matmul_precision('high')

    # # API NOT FINAL
    # # default: optimizes for large models, low compile-time
    # #          and no extra memory usage
    #compiled_model = torch.compile(model)
    #
    # reduce-overhead: optimizes to reduce the framework overhead
    #                and uses some extra memory. Helps speed up small models
    # compiled_model = torch.compile(model, mode="reduce-overhead")

    # # max-autotune: optimizes to produce the fastest model,
    # #               but takes a very long time to compile
    #compiled_model = torch.compile(model, mode="max-autotune")#加速

    compiled_model = model#model.half().float()

    best_model_loss = 50000
    best_model_both_auc = 0
    train_logprog = tqdm(range(epoch, args.epochs), dynamic_ncols=True)
    # torch.save(model.state_dict(), '/mnt/e/deep_learning/hyperspec_rgb_photo/ckpt/model_0.pt')
    scaler = torch.cuda.amp.GradScaler(enabled=args.amp == 1)
    #scaler = torch.cuda.amp.GradScaler()
    
    for epoch in train_logprog:
        model.train()
        outputs = []
        labels = []
        train_loss = []
        train_auc = []
        for step,data in enumerate(tqdm(train_loader)):
            # print(data)
            input = data['data'].to(device=device, dtype=torch.float32)
            input = input.to(torch.long)
            #ic(input,data_feature)
            label = data['Label'].to(device=device, dtype=torch.float32)
            if args.flag:
                label = label.to(device)
            else:
                label = label.to(device).view(-1, 1)
            #continue
            with torch.cuda.amp.autocast(dtype=torch.float16,enabled=args.amp==1):#混合精度
                # from torch.profiler import profile, ProfilerActivity, record_function
                # with profile(activities=[ProfilerActivity.CPU, ProfilerActivity.CUDA],record_shapes=True) as prof:
                #     with record_function("model_inference"):
                if args.model == 'Graphormer2Link_gate_type':
                    output = compiled_model(data_feature,train_adj,input,att_bias)
                    #ic(output)
                else:
                    x0,output = compiled_model(data_feature,train_adj,input,gene_num)
                #print(prof.key_averages().table(sort_by="cpu_time_total", row_limit=10))
                # np.savez(f'/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/graphormer2Link_result/Graphormer1Link_gate+/{step}_x0.npz',
                #         x0=(x0.cpu().detach().numpy()))
                if args.flag:
                    output = torch.softmax(output, dim=1)
                else:
                    output = torch.sigmoid(output)
                    output = output.squeeze(1)
                
                # loss_fn = nn.BCEWithLogitsLoss()
                # loss = loss_fn(input=output, target=label.float())
                #ic(output.shape)
                loss = F.binary_cross_entropy(output, label)
                    # 检查输出数据
                # criter = torch.nn.BCEWithLogitsLoss()
                # loss=criter(output, label)
                output=output.detach().cpu().float().numpy()
                label=label.cpu().float()
                #auc = roc_auc_score(label,output)
            #loss.backward(retain_graph=True)
            loss.backward()
            optimizer.step()
            #ic('=====backward==========')
            # scaler.scale(loss).backward()
            # scaler.step(optimizer)
            # scaler.update()

            optimizer.zero_grad(set_to_none=True)#set to none faster
            #scheduler.step()
            if args.lr_scheduler in ['CyclicLR','OneCycleLR','CosineAnnealingWarmRestarts']:
                #pr
                scheduler.step()
            if args.lr_scheduler not in ['CyclicLR', 'OneCycleLR', 'CosineAnnealingWarmRestarts']:
                if args.lr_scheduler in ['ReduceLROnPlateau']:
                    scheduler.step(metrics=np.mean(train_loss))
                else:
                    scheduler.step()
            outputs.append(output)
            labels.append(label)
            train_loss.append(np.array(loss.cpu().detach().float().numpy()))#loss

        outputs = np.concatenate(outputs)  # n,3
        labels = np.concatenate(labels)  # n,3
        #ic(labels,outputs)
        #ic(labels.shape,outputs.shape)
        train_avg_loss = np.mean(train_loss)
        train_avg_auc = roc_auc_score(labels,outputs)
        train_logprog.set_postfix(epoch=epoch, all_loss=train_avg_loss,all_auc=train_avg_auc)
        # print(outputs.shape,labels.shape)
        if (epoch + 1) % 50 == 0:
            avg_auc = roc_auc_score(labels, outputs)
            logger.info(
                f"train epoch:{epoch},avg_auc,{avg_auc}")

        val_outputs,val_labels,val_score_dict = evaluation(args, model, val_loader,data_feature,train_adj,gene_num,att_bias)

        # val_outputs = [i for i in val_outputs[:,0]]
        # ic(val_outputs,val_labels)
        # val_loss = [F.binary_cross_entropy(val_outputs[i], val_labels[i]) for i in range(len(val_outputs))]
        # val_loss_mean = np.mean(val_loss)
        # tb.add_scalar("val_loss", val_loss_mean, epoch)

        test_outputs,test_labels,test_score_dict = evaluation(args, model, test_loader,data_feature,train_adj,gene_num,att_bias) 
 
        # test_loss = mean_squared_error(test_outputs, test_labels) #测试集损失随epoch变化
        # test_loss_mean = np.mean(test_loss) #测试集损失随epoch变化
        

        avg_both_auc = (val_score_dict['AUC'] + test_score_dict['AUC']) / 2  # val_r2

        if avg_both_auc > best_model_both_auc:
            # print(outputs[:,0]-labels[:,0])
            best_model_both_auc = avg_both_auc
            torch.save(
                {"state": model.state_dict(),  # 模型参数
                 "epoch": epoch
                 }
                , 'best_model_both.pt')
            # logger.info(f"val epoch:{epoch},best_model_r2,{best_model_r2}")
            logger.info(f"val epoch:{epoch},val_auc,{val_score_dict['AUC']}")  # ,
            logger.info(f"test_auc:{test_score_dict['AUC']}")


    logger.info(f"---------------------------------best model results------------------------------------------")  #
    model.load_state_dict(torch.load('best_model_both.pt', map_location=device)['state'])
    train_outputs,train_labels,train_score_dict = evaluation(args, model, train_loader,data_feature,train_adj,gene_num,att_bias)
    np.savez('train_out.npz',
             labels=(train_labels),  #
             preds=(train_outputs),  #
             )
    logger.info(train_score_dict)
        #f"train_AUC: {train_score_dict['AUC']:.3f} train_AUPR: {train_score_dict['AUPR']:.3f} train_ACC: {train_score_dict['ACC']:.3f} train_F1_score: {train_score_dict['F1_score']:.3f} train_brier_score: {train_score_dict['brier_score']:.3f} train_mr: {train_score_dict['MR']:.3f} train_mrr: {train_score_dict['MRR']:.3f} train_hits_at_1: {train_score_dict['hits@1']:.3f} train_hits_at_3: {train_score_dict['hits@3']:.3f} train_hits_at_10: {train_score_dict['hits@10']:.3f}")
    val_outputs,val_labels,val_score_dict = evaluation(args, model, val_loader,data_feature,train_adj,gene_num,att_bias)
    np.savez('val_out.npz',
             labels=(val_labels),  #
             preds=(val_outputs),  #
             )
    logger.info(val_score_dict)
        #f"val_AUC: {val_AUC:.3f} val_AUPR: {val_AUPR:.3f} val_ACC: {val_ACC:.3f} val_F1_score: {val_F1_score:.3f} val_brier_score: {val_brier_score:.3f}  val_mr: {val_score_dict['MR']:.3f} val_mrr: {val_score_dict['MRR']:.3f} val_hits_at_1: {val_score_dict['hits@1']:.3f} val_hits_at_3: {val_score_dict['hits@3']:.3f} val_hits_at_10: {val_score_dict['hits@10']:.3f}")

    test_outputs,test_labels,test_score_dict = evaluation(args, model, test_loader,data_feature,train_adj,gene_num,att_bias)
    np.savez('test_out.npz',
             labels=(test_labels),  #
             preds=(test_outputs),  #
             )
    logger.info(test_score_dict)
        #f"test_AUC: {test_AUC:.3f} test_AUPR: {test_AUPR:.3f} test_ACC: {test_ACC:.3f} test_F1_score: {test_F1_score:.3f} test_brier_score: {test_brier_score:.3f} test_mr: {test_score_dict['MR']:.3f} test_mrr: {test_score_dict['MRR']:.3f} test_hits_at_1: {test_score_dict['hits@1']:.3f} test_hits_at_3: {test_score_dict['hits@3']:.3f} test_hits_at_10: {test_score_dict['hits@10']:.3f}")

    result = [args.name,(val_score_dict['AUC']+test_score_dict['AUC'])/2,(val_score_dict['AUPR']+test_score_dict['AUPR'])/2,(val_score_dict['ACC']+test_score_dict['ACC'])/2,(val_score_dict['F1_score']+test_score_dict['F1_score'])/2,(val_score_dict['brier_score']+test_score_dict['brier_score'])/2]
    result = pd.DataFrame(result).transpose()
    result.to_csv(args.expression.output_path, mode='a', header=False, index=False)
    # with open('/mnt/h/study/deep_learning/outputs/gene_result.csv', mode='a') as f:
    #     f.write(
    #         f'{args.name},{(val_score_dict['AUC']+test_score_dict['AUC'])/2},{(val_score_dict['AUPR']+test_score_dict['AUPR'])/2},{(val_score_dict['ACC']+test_score_dict['ACC'])/2},{(val_score_dict['F1_score']+test_score_dict['F1_score'])/2},{(val_score_dict['brier_score']+test_score_dict['brier_score'])/2} \n')
            
            
            
            

    # #
    # logger.info(f"lr:{lr},weight_decay:{weight_decay}")  #
    # history = {
    #     "lr": lr,
    #     "weight_decay": weight_decay,
    #     "train_Vcmax_r2": train_Vcmax_r2,
    #     "train_Jmax_r2": train_Jmax_r2,
    #     "val_Vcmax_r2": val_Vcmax_r2,
    #     "val_Jmax_r2": val_Jmax_r2,
    #     "test_Vcmax_r2": test_Vcmax_r2,
    #     "test_Jmax_r2": test_Jmax_r2,
    #     "optim": args.optim,
    #     "lr_scheduler": args.lr_scheduler,
    #     "batch_size": args.batch_size
    # }
    # json.dump(history, open(args.history_file, "a"), indent=2)


# @hydra.main(config_path="conf", config_name='config.yaml')
# def main(args):
#     try:
#         # print(args)
#         _main(args)

#     except Exception as e:
#         logger.exception(e)
#         # Hydra intercepts exit code, fixed in beta but I could not get the beta to work
#         os._exit(1)


from sacred import Experiment
from sacred.observers import MongoObserver
from sacred.stflow import LogFileWriter
import sacred

from omegaconf import DictConfig
@hydra.main(config_path="conf", config_name='config.yaml')
def main(cfg:DictConfig):
    #try:
        # if cfg.sacred:
        #     ex = Experiment(cfg.name)
        #     ex.observers.append(MongoObserver.create(url='localhost:27017',
        #                                              db_name='sacred'))
        #     ex.add_config({'_args':cfg})
        #     @ex.main
        #     #@LogFileWriter(ex)
        #     def run(_args,_run):
        #         return _main(_args,_run)

        #     ex.run()
        # else:
    _main(cfg)
    # except Exception as e:
    #     #logger.exception("Some error happened")
    #     torch.cuda.empty_cache()
    #     logger.exception(e)
    #     # Hydra intercepts exit code, fixed in beta but I could not get the beta to work
    #     os._exit(1)



if __name__ == "__main__":
    # import os
    # os.environ['CUDA_LAUNCH_BLOCKING'] = '1'
    main()
