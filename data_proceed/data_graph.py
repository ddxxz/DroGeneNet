import torch
import numpy as np
from torch.utils.data import Dataset
from scipy import signal
from torch.utils.data import DataLoader
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler
import scipy.sparse as sp
import numpy as np
from sklearn.metrics import roc_auc_score,average_precision_score
import torch.nn as nn
from sklearn import metrics
from icecream import install , ic
install()
import random
class BasicDataset(Dataset):
    def __init__(self,args,train_set,num_gene,flag=False):
        self.train_set = train_set
        self.train_data = pd.read_csv(train_set, index_col=0).values
        self.num_gene = num_gene
        self.flag = flag
        self.train_input = self.train_data[:,:2]
        self.train_label = self.train_data[:,-1]
        #ic(len(train_label))
        if self.flag:
            train_len = len(self.train_label)           # a  <-> b
            train_tan = np.zeros([train_len,2])
            train_tan[:,0] = 1 - self.train_label       #0  a - b    b -> a 
            train_tan[:,1] = self.train_label           #1  a -> b   b - a
            self.train_label = train_tan  #转化为单向

    def __len__(self):
        return len(self.train_data)

    def __getitem__(self, idx):
        

        data = self.train_input[idx].astype(np.int64)
        label = self.train_label[idx].astype(np.float32)

        dataset = {'data': torch.from_numpy(np.array(data, dtype='int')),
                    'Label': torch.from_numpy(np.array(label, dtype='float32'))
                    }
        return dataset
    def Adj_Generate(self,TF_set,adj_type,flag=False, loop=False):
        adj = sp.dok_matrix((self.num_gene, self.num_gene), dtype=np.float32)
        for pos in self.train_data:
            #ic(pos)
            tf = pos[0]
            target = pos[1]

            if flag == False:
                if pos[-1] == 1:
                    #ic(adj.shape)
                    adj[tf, target] = 1.0
                    adj[target, tf] = 1.0
            else:
                if pos[-1] == 1:
                    #ic(tf,target)
                    adj[tf, target] = 1.0
                    if target in TF_set:
                        adj[target, tf] = 1.0

        if loop:
            adj = adj + sp.identity(self.num_gene)

        if adj_type=='adj_drop':
            # 添加噪声
            keys = list(adj.keys())
            # 随机删除50%的边
            for key in keys:
                if random.random() < 0.5:  # 有50%的概率删除边
                    adj[key] = 0
        elif adj_type=='adj_weight':
            # 为存在的边赋予0.6的权重
            keys = list(adj.keys())
            for key in keys:
                if adj[key] != 0:
                    adj[key] = 0.6
            num_possible_edges = self.num_gene * (self.num_gene - 1) - len(adj)
            num_edges_to_select = int(0.00001 * num_possible_edges)
            ic('start')
            # 随机选择边并赋权重，不再需要提前列出所有可能的边
            ic(num_edges_to_select)
            for _ in range(num_edges_to_select):
                while True:
                    i, j = np.random.randint(self.num_gene), np.random.randint(self.num_gene)
                    if i != j and (i, j) not in adj:
                        adj[i, j] = 0.4
                        break

        adj = adj.todok()
        ic(adj.shape,len(adj.keys()))

        return adj

class load_data():
    def __init__(self, data, normalize=True):
        #data = pd.concat(data.values(), axis=1).to_numpy()
        self.data = data
        self.normalize = normalize

    def data_normalize(self,data):
        #data = pd.concat(data.values(), axis=1).to_numpy()
        standard = StandardScaler()
        epr = standard.fit_transform(data.T)
        return epr.T


    def exp_data(self):
        data_feature = self.data.values
        if self.normalize:
            data_feature = self.data_normalize(data_feature)
        data_feature = data_feature.astype(np.float32)
        return data_feature


def adj2saprse_tensor(adj):
    coo = adj.tocoo()
    i = torch.LongTensor(np.array([coo.row, coo.col]))
    v = torch.from_numpy(coo.data).float()
    ic(i.shape, i[:,0],v.shape,v[0])
    adj_sp_tensor = torch.sparse_coo_tensor(i, v, coo.shape)
    return adj_sp_tensor

def Evaluation(y_true, y_pred,flag=False):
    if flag:
        # y_p = torch.argmax(y_pred,dim=1)
        y_p = y_pred[:,-1]
        y_p = y_p.cpu().detach().numpy()
        y_p = y_p.flatten()
    else:
        y_p = y_pred.cpu().detach().numpy()
        y_p = y_p.flatten()

    y_pred_binary = [1 if p >= 0.5 else 0 for p in y_p]
    y_t = y_true.cpu().numpy().flatten().astype(int)

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
    log_loss = calculate_log_loss_weighted(y_t, y_p)
    log_loss = np.mean(log_loss)
    brier_score = calculate_brier_score(y_t, y_p)
    return AUC, AUPR, ACC,TP,TN,FP,FN,precision,recall,F1_score,log_loss,brier_score

def calculate_TP(y, y_pred): 
    tp = 0 
    for i, j in zip(y, y_pred): 
        if i == j == 1: 
            tp += 1 
    return tp 
def calculate_TN(y, y_pred): 
    tn = 0 
    for i, j in zip(y, y_pred): 
        if i == j == 0: 
            tn += 1 
    return tn 
def calculate_FP(y, y_pred): 
    fp = 0 
    for i, j in zip(y, y_pred): 
        if i == 0 and j == 1: 
            fp += 1 
    return fp 
def calculate_FN(y, y_pred): 
    fn = 0 
    for i, j in zip(y, y_pred): 
        if i == 1 and j == 0: 
            fn += 1 
    return fn
 
def calculate_precision(y, y_pred): 
    tp = calculate_TP(y, y_pred) 
    fp = calculate_FP(y, y_pred) 
    return tp / (tp + fp+1e-10)

def calculate_recall(y, y_pred): 
    tp = calculate_TP(y, y_pred) 
    fn = calculate_FN(y, y_pred) 
    return tp / (tp + fn+1e-10)

def calculate_F1(y, y_pred): 
    p = calculate_precision(y, y_pred) 
    r = calculate_recall(y, y_pred) 
    return 2*p*r / (p+r+1e-10)

def calculate_log_loss(t, p): 
    log_loss = -1.0*(t*np.log(p) + (1-t)*(t*np.log(1-p)))
    return log_loss

def calculate_log_loss_weighted(t, p): 
    w1=0.33
    w2=0.67
    log_loss = -1.0*(w1*t*np.log(p) + w2*(1-t)*(t*np.log(1-p))) 
    return log_loss

def calculate_brier_score(y, y_pred): 
    s=0 
    for i, j in zip(y, y_pred): 
        s += (j-i)**2 
    return s * (1/len(y))

def calculate_MRR_HITS(y, y_pred):
    # 对预测的分数降序排列，得到索引
    sorted_idx = np.argsort(y_pred)[::-1]

    # 将实际标签按照对应的排序索引进行重新排序
    sorted_true = y[sorted_idx]

    # 计算每个正例的排名（即第一个出现的位置，加1因为排名从1开始计数）
    ranks = np.where(sorted_true == 1)[0] + 1
    # 计算 MR
    mr = np.argwhere(sorted_true == 1)[0][0] + 1
    # 计算 MRR
    mrr = 1 / (np.argwhere(sorted_true == 1)[0][0] + 1)
    # 计算 HITS@1 和 HITS@3
    # 如果排名在top N内则为命中
    hits_at_1 = np.mean(ranks <= 1)
    hits_at_3 = np.mean(ranks <= 3)
    hits_at_10 = np.mean(ranks <= 10)

    return mr,mrr,hits_at_1,hits_at_3,hits_at_10

def safe_ranking(_scores, verbose=False):
    import random
    import numpy as np

    if _scores is None:  # 检查_scores是否为None
        raise ValueError("Input array _scores is empty")
    #ic(_scores)
    pos_score = _scores[0]
    same_score_loc = np.where(_scores == pos_score)[0]
    assert same_score_loc.size > 0 and same_score_loc[0] == 0
    rdm_pos_loc = same_score_loc[random.randint(0, same_score_loc.shape[0] - 1)]
    _sort_idxs = np.argsort(-_scores)
    _rank = np.where(_sort_idxs == rdm_pos_loc)[0][0] + 1

    if verbose:
        _default_rank = np.where(_sort_idxs == 0)[0][0] + 1
        print("From safe_ranking: default rank is {}, after safe {}".format(_default_rank, _rank))

    return _rank

def calculate_MR_MRR_HITS(y, y_pred):
    scores = y_pred
    ranks = np.zeros(len(scores))  # 用于存储每个样本的排名
    hits = np.zeros((1000, len(scores)))  # 用于存储hits指标
    for i in range(len(scores)):
        ranks[i] = safe_ranking(scores)  # 计算每个样本的排名
        for hits_level in range(1000):
            if ranks[i] <= hits_level + 1:
                hits[hits_level, i] = 1.0
            else:
                hits[hits_level, i] = 0.0
    
    score_dict = {
        'hits@1': np.mean(hits[0]),
        'hits@3': np.mean(hits[2]),
        'hits@10': np.mean(hits[9]),
        'hits@30': np.mean(hits[29]),
        'hits@50': np.mean(hits[49]),
        'hits@100': np.mean(hits[99]),
        'hits@300': np.mean(hits[299]),
        'hits@500': np.mean(hits[499]),
        'hits@1000': np.mean(hits[999]),
        'MR': np.mean(ranks),
        'MRR': np.mean(1. / ranks),
        'Q95R': np.quantile(ranks, q=0.95),
    }

    return score_dict

def calculate_hits(scores, k):
    if isinstance(scores, list):
        sorted_indices = sorted(range(len(scores)), key=lambda x: scores[x], reverse=True)
    elif isinstance(scores, dict):
        sorted_indices = [key for key, value in sorted(scores.items(), key=lambda x: x[1]['prediction'], reverse=True)]
    elif isinstance(scores, np.ndarray):
        sorted_indices = np.argsort(-scores)
    elif isinstance(scores, pd.Series):
        sorted_indices = scores.sort_values(ascending=False).index.tolist()
    else:
        raise ValueError("Unsupported scores format")

    hits_at_k = sum(1 for idx in sorted_indices[:k] if scores[idx] >= 1.0) / k

    return hits_at_k
def calculate_mr_mrr(scores):
    if isinstance(scores, list):
        sorted_indices = sorted(range(len(scores)), key=lambda x: scores[x])
    elif isinstance(scores, dict):
        sorted_indices = [key for key, value in sorted(scores.items(), key=lambda x: x[1]['prediction'])]
    elif isinstance(scores, np.ndarray):
        sorted_indices = np.argsort(scores)
    elif isinstance(scores, pd.Series):
        sorted_indices = scores.sort_values().index.tolist()
    else:
        raise ValueError("Unsupported scores format")

    ranks = np.arange(1, len(sorted_indices) + 1)
    mr = np.mean(ranks)
    mrr = np.mean(1.0 / ranks)

    return mr, mrr

def normalize(expression):
    std = StandardScaler()
    epr = std.fit_transform(expression)

    return epr

def Network_Statistic(data_type,net_scale,net_type):

    if net_type =='STRING':
        dic = {'hESC500': 0.024, 'hESC1000': 0.021, 'hHEP500': 0.028, 'hHEP1000': 0.024, 'mDC500': 0.038,
               'mDC1000': 0.032, 'mESC500': 0.024, 'mESC1000': 0.021, 'mHSC-E500': 0.029, 'mHSC-E1000': 0.027,
               'mHSC-GM500': 0.040, 'mHSC-GM1000': 0.037, 'mHSC-L500': 0.048, 'mHSC-L1000': 0.045}

        query = data_type + str(net_scale)
        scale = dic[query]
        return scale



    elif net_type == 'Non-Specific':

        dic = {'hESC500': 0.016, 'hESC1000': 0.014, 'hHEP500': 0.015, 'hHEP1000': 0.013, 'mDC500': 0.019,
               'mDC1000': 0.016, 'mESC500': 0.015, 'mESC1000': 0.013, 'mHSC-E500': 0.022, 'mHSC-E1000': 0.020,
               'mHSC-GM500': 0.030, 'mHSC-GM1000': 0.029, 'mHSC-L500': 0.048, 'mHSC-L1000': 0.043}

        query = data_type + str(net_scale)
        scale = dic[query]
        return scale

    elif net_type == 'Specific':
        dic = {'hESC500': 0.164, 'hESC1000': 0.165,'hHEP500': 0.379, 'hHEP1000': 0.377,'mDC500': 0.085,
               'mDC1000': 0.082,'mESC500': 0.345, 'mESC1000': 0.347,'mHSC-E500': 0.578, 'mHSC-E1000': 0.566,
               'mHSC-GM500': 0.543, 'mHSC-GM1000': 0.565,'mHSC-L500': 0.525, 'mHSC-L1000': 0.507}

        query = data_type + str(net_scale)
        scale = dic[query]
        return scale

    elif net_type == 'Lofgof':
        dic = {'mESC500': 0.158, 'mESC1000': 0.154}

        query = 'mESC' + str(net_scale)
        scale = dic[query]
        return scale

    else:
        raise ValueError
    
if __name__ == '__main__':
    data = BasicDataset('data/rgb_read_noseg.csv')#实例化
    dataloader = DataLoader(data,batch_size=1)
    #for data in dataloader:
    generator = iter(dataloader)
    data = next(generator)
    spectral = data['Spectral']
    rgb = data['RGB']
    label = data['Chl']
    print(spectral.shape)
    print(rgb.shape)
    print(label.shape)
    print(label)
        #break
    data = next(generator)

