import torch
import torch.nn as nn
import numpy as np
import sys
import pandas as pd
import torch.nn.functional as F
import torch.optim as optm
from torch.nn import CosineSimilarity
from icecream import install , ic
from torch.nn import ModuleList, Linear
from icecream import install , ic
install()
sys.path.append('/home/win/4T/GeneDL/OSDrought_GCNAT_Link')
from models.BaseModel import BaseModel
import torch.nn.functional as F


class GENELink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, 
                  flag, device, type, hidden_dim,visualize,**kwgs):
        super(GENELink, self).__init__(device)
        self.num_head1 = 3
        self.num_head2 = 3
        self.device = device
        self.alpha = 0.2
        self.type = type
        self.flag = flag
        self.reduction = 'concat'

        if self.reduction == 'mean':
            self.hidden1_dim = hidden1_dim
            self.hidden2_dim = hidden2_dim
        elif self.reduction == 'concat':
            self.hidden1_dim = self.num_head1*hidden1_dim
            self.hidden2_dim = self.num_head2*hidden2_dim

        # self.layer_options =nn.ModuleDict({
        #     'japonica_CK': nn.Linear(602, 1000).to(device),
        #     'japonica_CD': nn.Linear(635, 1000).to(device),
        #     'japonica_CE': nn.Linear(85, 1000).to(device),
        #     'indica_CK': nn.Linear(333, 1000).to(device),
        #     'indica_CD': nn.Linear(263, 1000).to(device),
        #     'indica_CE': nn.Linear(56, 1000).to(device),
        #     'zeamays_CK': nn.Linear(743, 1000).to(device),
        #     'zeamays_CD': nn.Linear(777, 1000).to(device),
        #     'zeamays_CE': nn.Linear(10, 1000).to(device),
        #     'homo_CK': nn.Linear(1676, 1000).to(device),
        #     'homo_CD': nn.Linear(1675, 1000).to(device),
        #     'homo_CE': nn.Linear(151, 1000).to(device),
        # })

        self.ConvLayer1 = [AttentionLayer(input_dim,hidden1_dim,self.alpha) for _ in range(self.num_head1)]

        for i, attention in enumerate(self.ConvLayer1):
            self.add_module('ConvLayer1_AttentionHead{}'.format(i),attention)

        self.ConvLayer2 = [AttentionLayer(self.hidden1_dim,hidden2_dim,self.alpha) for _ in range(self.num_head2)]

        for i, attention in enumerate(self.ConvLayer2):
            self.add_module('ConvLayer2_AttentionHead{}'.format(i),attention)

        self.tf_linear1 = nn.Linear(hidden2_dim,hidden3_dim)
        self.target_linear1 = nn.Linear(hidden2_dim,hidden3_dim)

        self.tf_linear2 = nn.Linear(hidden3_dim,output_dim)
        self.target_linear2 = nn.Linear(hidden3_dim, output_dim)

        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Linear(2*output_dim, 2)

        self.reset_parameters()

    def reset_parameters(self):
        for attention in self.ConvLayer1:
            attention.reset_parameters()

        for attention in self.ConvLayer2:
            attention.reset_parameters()

        nn.init.xavier_uniform_(self.tf_linear1.weight,gain=1.414)
        nn.init.xavier_uniform_(self.target_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.tf_linear2.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear2.weight, gain=1.414)

    def encode(self,x,adj):
        # ic(x.shape,adj.shape)
        # #---------适合有向图，取所有边----------------------------
        # if self.flag:
        #     adj_matrix_coalesced = adj.coalesce()

        #     # 获取稀疏张量的坐标
        #     sparse_coords = adj_matrix_coalesced.indices()

        #     # 有向图的边不需要上三角筛选
        #     edge_index = sparse_coords
        # #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        # else:
        #     adj_matrix=adj
        #     adj_matrix_coalesced = adj_matrix.coalesce()
        #     #ic(adj_matrix_coalesced)
        #     # 获取聚合稀疏张量的坐标和值
        #     sparse_coords = adj_matrix_coalesced.indices()
        #     values = adj_matrix_coalesced.values()
        #     # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
        #     row_indices, col_indices = sparse_coords[0], sparse_coords[1]
        #     mask = row_indices < col_indices
        #     edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        # ic(edge_index.shape)

        if self.reduction =='concat':
            x = torch.cat([att(x, adj) for att in self.ConvLayer1], dim=1)
            x = F.elu(x)
        elif self.reduction =='mean':
            x = torch.mean(torch.stack([att(x, adj) for att in self.ConvLayer1]), dim=0)
            x = F.elu(x)
        else:
            raise TypeError

        x = torch.mean(torch.stack([att(x, adj) for att in self.ConvLayer2]),dim=0)

        return x  #50,32

    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))


    def forward(self,x,adj,train_sample,gene_num):
        #ic(x.shape,adj.shape,train_sample.shape)#x 910,758  adj 910,910  train_sample 256 2
        #ic(x.shape)
        # if x.shape[1] == 602:
        #     x = self.layer_options['japonica_CK'](x)
        # elif x.shape[1] == 635:
        #     x = self.layer_options['japonica_CD'](x)
        # elif x.shape[1] == 85:
        #     x = self.layer_options['japonica_CE'](x)
        # elif x.shape[1] == 333:
        #     x = self.layer_options['indica_CK'](x)
        # elif x.shape[1] == 263:
        #     x = self.layer_options['indica_CD'](x)
        # elif x.shape[1] == 56:
        #     x = self.layer_options['indica_CE'](x)
        # elif x.shape[1] == 743:
        #     x = self.layer_options['zeamays_CK'](x)
        # elif x.shape[1] == 777:
        #    # print(self.layer_options['zeamays_CD'])
        #     x = self.layer_options['zeamays_CD'](x)
        # elif x.shape[1] == 10:
        #     x = self.layer_options['zeamays_CE'](x)

        # elif x.shape[1] == 1676:
        #     x = self.layer_options['homo_CK'](x)
        # elif x.shape[1] == 1675:
        #    # print(self.layer_options['zeamays_CD'])
        #     x = self.layer_options['homo_CD'](x)
        # elif x.shape[1] == 151:
        #     x = self.layer_options['homo_CE'](x)

        embed = self.encode(x,adj) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
        train_tf = tf_embed[train_sample[:,0]]
        #ic(train_tf.shape)#train_tf 256 16
        train_target = target_embed[train_sample[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output

class AttentionLayer(BaseModel):
    def __init__(self,input_dim,output_dim,alpha=0.2,bias=True):
        super(AttentionLayer, self).__init__()
        self.input_dim = input_dim
        self.output_dim = output_dim
        self.alpha = alpha
        self.weight = nn.Parameter(torch.FloatTensor(self.input_dim, self.output_dim))
        self.weight_interact = nn.Parameter(torch.FloatTensor(self.input_dim,self.output_dim))
        self.a = nn.Parameter(torch.zeros(size=(2*self.output_dim,1)))
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(self.output_dim))
        else:
            self.register_parameter('bias', None)

        self.reset_parameters()
    def reset_parameters(self):
        nn.init.xavier_uniform_(self.weight.data, gain=1.414)
        nn.init.xavier_uniform_(self.weight_interact.data, gain=1.414)
        if self.bias is not None:
            self.bias.data.fill_(0)
        nn.init.xavier_uniform_(self.a.data, gain=1.414)

    def _prepare_attentional_mechanism_input(self, x):

        Wh1 = torch.matmul(x, self.a[:self.output_dim, :])
        Wh2 = torch.matmul(x, self.a[self.output_dim:, :])
        e = F.leaky_relu(Wh1 + Wh2.T,negative_slope=self.alpha)
        return e

    def forward(self,x,adj):

        #print(x.shape)
        h = torch.matmul(x, self.weight)
        #ic(h.shape)
        e = self._prepare_attentional_mechanism_input(h)
        #ic(e.shape)
        zero_vec = -9e15 * torch.ones_like(e)
        attention = torch.where(adj.to_dense()>0, e, zero_vec)
        attention = F.softmax(attention, dim=1)
        # attention = F.softmax(e, dim=1)

        attention = F.dropout(attention, training=self.training)
        h_pass = torch.matmul(attention, h)

        output_data = h_pass
        output_data = F.leaky_relu(output_data,negative_slope=self.alpha)
        output_data = F.normalize(output_data,p=2,dim=1)
        if self.bias is not None:
            output_data = output_data + self.bias
        return output_data

class GraphConvolution(BaseModel):
    """
    简单的GCN层，相似于 https://arxiv.org/abs/1609.02907
    """
    def __init__(self, in_features, out_features, bias=True):
        super(GraphConvolution, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = nn.Parameter(torch.FloatTensor(in_features, out_features))
        if bias:
            self.bias = nn.Parameter(torch.FloatTensor(out_features))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        nn.init.kaiming_uniform_(self.weight)
        if self.bias is not None:
            nn.init.zeros_(self.bias)

    def forward(self, input, adj):
        #ic(input.shape)
        support = torch.mm(input, self.weight)
        output = torch.mm(adj, support)
        #ic(output.shape)
        if self.bias is not None:
            return output + self.bias
        else:
            return output

class GNNLink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, 
                  flag, device, type, hidden_dim,visualize,**kwgs):
        super(GNNLink, self).__init__(device)
        self.device = device
        self.num_head1 = 3
        self.num_head2 = 3
        self.reduction = 'concat'
        self.type = type
        self.flag = flag
        # 对于GCN，我们不使用多头机制
        self.gcn1 = GraphConvolution(input_dim, hidden1_dim)
        self.gcn2 = GraphConvolution(hidden1_dim, hidden1_dim)

        # 对于GAT，我们根据头数创建多个注意力机制层
        # self.attention_heads = ModuleList([
        #     AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        # ])

        # self.attention_heads = ModuleList([
        #     AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        # ])
        # 如果采用拼接，更新hidden2_dim值为所有头输出拼接后的维度
        if self.reduction == 'concat':
            hidden2_dim = self.num_head2 * (hidden2_dim)
        elif self.reduction != 'mean':
            raise ValueError("Reduction method must be 'mean' or 'concat'.")

        self.hidden2_dim = hidden2_dim
        
        self.tf_linear1 = Linear(128, hidden3_dim)
        self.target_linear1 = Linear(128, hidden3_dim)

        self.tf_linear2 = Linear(hidden3_dim, output_dim)
        self.target_linear2 = Linear(hidden3_dim, output_dim)

        # 我们可以增加一个分类头，如果我们执行一个节点对分类任务
        self.linear = Linear(2 * output_dim, 2)
        self.decode_linear = Linear(output_dim, 2)
        
        self.reset_parameters()

    def reset_parameters(self):
        self.gcn1.reset_parameters()
        self.gcn2.reset_parameters()
        # for head in self.attention_heads:
        #     head.reset_parameters()

        # 初始化其余的全连接层
        nn.init.xavier_uniform_(self.tf_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.tf_linear2.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear2.weight, gain=1.414)
        
    def encode(self, x, adj):
        # #---------适合有向图，取所有边----------------------------
        # if self.flag:
        #     adj_matrix_coalesced = adj.coalesce()

        #     # 获取稀疏张量的坐标
        #     sparse_coords = adj_matrix_coalesced.indices()

        #     # 有向图的边不需要上三角筛选
        #     edge_index = sparse_coords
        # #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        # else:
        #     adj_matrix=adj
        #     adj_matrix_coalesced = adj_matrix.coalesce()
        #     #ic(adj_matrix_coalesced)
        #     # 获取聚合稀疏张量的坐标和值
        #     sparse_coords = adj_matrix_coalesced.indices()
        #     values = adj_matrix_coalesced.values()
        #     # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
        #     row_indices, col_indices = sparse_coords[0], sparse_coords[1]
        #     mask = row_indices < col_indices
        #     edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        # adj=edge_index
        # GCN 编码
        x = F.relu(self.gcn1(x, adj))
        x = F.relu(self.gcn2(x, adj))
        #ic('gcn',x.shape)
        # GAT 编码
        #x_all = [att(x, adj) for att in self.attention_heads]
        #ic(x_all[1].shape)
        # if self.reduction == 'concat':
        #     x = torch.cat(x, dim=1)
        # elif self.reduction == 'mean':
        #     x = torch.mean(torch.stack(x, dim=0), dim=0)
        # else:
        #     raise ValueError(f"Unknown reduction: {self.reduction}")
        #ic('out',x.shape)
        return x  #50,32
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
        
    def forward(self,x,adj,train_sample,gene_num):
        #ic(x.shape,adj.shape,train_sample.shape)#x 910,758  adj 910,910  train_sample 256 2
        #ic('-----',x.shape)
        embed = self.encode(x,adj) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
        train_tf = tf_embed[train_sample[:,0]]
        #ic(train_tf.shape)#train_tf 256 16
        train_target = target_embed[train_sample[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output


class GCNGATLink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, num_head1, num_head2,
                 alpha, device, type,reduction='mean'):
        super(GCNGATLink, self).__init__(device)
        self.device = device
        self.num_head1 = num_head1
        self.num_head2 = num_head2
        self.reduction = reduction
        self.type = type
        
        # 对于GCN，我们不使用多头机制
        self.gcn1 = GraphConvolution(input_dim, hidden1_dim)
        #self.gcn2 = GraphConvolution(hidden1_dim, hidden1_dim)

        # 对于GAT，我们根据头数创建多个注意力机制层
        # self.attention_heads = ModuleList([
        #     AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        # ])

        self.attention_heads = ModuleList([
            AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        ])
        # 如果采用拼接，更新hidden2_dim值为所有头输出拼接后的维度
        if reduction == 'concat':
            hidden2_dim = num_head2 * (hidden2_dim)
        elif reduction != 'mean':
            raise ValueError("Reduction method must be 'mean' or 'concat'.")

        self.hidden2_dim = hidden2_dim
        
        self.tf_linear1 = Linear(hidden2_dim, hidden3_dim)
        self.target_linear1 = Linear(hidden2_dim, hidden3_dim)

        self.tf_linear2 = Linear(hidden3_dim, output_dim)
        self.target_linear2 = Linear(hidden3_dim, output_dim)

        # 我们可以增加一个分类头，如果我们执行一个节点对分类任务
        self.linear = Linear(2 * output_dim, 2)
        
        self.reset_parameters()

    def reset_parameters(self):
        self.gcn1.reset_parameters()
        #self.gcn2.reset_parameters()
        for head in self.attention_heads:
            head.reset_parameters()

        # 初始化其余的全连接层
        nn.init.xavier_uniform_(self.tf_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.tf_linear2.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear2.weight, gain=1.414)
        
    def encode(self, x, adj):
        # GCN 编码
        x = F.relu(self.gcn1(x, adj))
        #x = F.relu(self.gcn2(x, adj))
        #ic('gcn',x.shape)
        # GAT 编码
        x_all = [att(x, adj) for att in self.attention_heads]
        #ic(x_all[1].shape)
        if self.reduction == 'concat':
            x = torch.cat(x_all, dim=1)
        elif self.reduction == 'mean':
            x = torch.mean(torch.stack(x_all, dim=0), dim=0)
        else:
            raise ValueError(f"Unknown reduction: {self.reduction}")
        #ic('out',x.shape)
        return x  #50,32
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            prob = torch.sum(prob,dim=1).view(-1,1)
            return prob
        elif self.type =='cosine':
            prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            h = torch.cat([tf_embed, target_embed],dim=1)
            prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
        
    def forward(self,x,adj,train_sample):
        #ic(x.shape,adj.shape,train_sample.shape)#x 910,758  adj 910,910  train_sample 256 2
        #ic('-----',x.shape)
        embed = self.encode(x,adj) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
        train_tf = tf_embed[train_sample[:,0]]
        #ic(train_tf.shape)#train_tf 256 16
        train_target = target_embed[train_sample[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output
            
class LargeGCNGATLink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, num_head1, num_head2,
                 alpha, device, type,reduction='mean'):
        super(LargeGCNGATLink, self).__init__(device)
        self.device = device
        self.num_head1 = num_head1
        self.num_head2 = num_head2
        self.reduction = reduction
        self.type = type
        
        # 对于GCN，我们不使用多头机制
        self.gcn1 = GraphConvolution(input_dim, hidden1_dim)
        #self.gcn2 = GraphConvolution(hidden1_dim, hidden1_dim)

        # 对于GAT，我们根据头数创建多个注意力机制层
        # self.attention_heads = ModuleList([
        #     AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        # ])

        self.attention_heads = ModuleList([
            AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        ])
        # 如果采用拼接，更新hidden2_dim值为所有头输出拼接后的维度
        if reduction == 'concat':
            hidden2_dim = num_head2 * (hidden2_dim)
        elif reduction != 'mean':
            raise ValueError("Reduction method must be 'mean' or 'concat'.")

        self.hidden2_dim = hidden2_dim
        
        self.tf_linear1 = Linear(hidden2_dim, hidden3_dim)
        self.target_linear1 = Linear(hidden2_dim, hidden3_dim)

        self.tf_linear2 = Linear(hidden3_dim, output_dim)
        self.target_linear2 = Linear(hidden3_dim, output_dim)

        # 我们可以增加一个分类头，如果我们执行一个节点对分类任务
        self.linear = Linear(2 * output_dim, 2)
        
        self.reset_parameters()

    def reset_parameters(self):
        self.gcn1.reset_parameters()
        #self.gcn2.reset_parameters()
        for head in self.attention_heads:
            head.reset_parameters()

        # 初始化其余的全连接层
        nn.init.xavier_uniform_(self.tf_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.tf_linear2.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear2.weight, gain=1.414)
        
    def encode(self, x, adj):
        # GCN 编码

        x = F.relu(self.gcn1(x, adj))
        #x = F.relu(self.gcn2(x, adj))
        #ic('gcn',x.shape)
        # GAT 编码
        x_all = [att(x, adj) for att in self.attention_heads]
        #ic(x_all[1].shape)
        if self.reduction == 'concat':
            x = torch.cat(x_all, dim=1)
        elif self.reduction == 'mean':
            x = torch.mean(torch.stack(x_all, dim=0), dim=0)
        else:
            raise ValueError(f"Unknown reduction: {self.reduction}")
        #ic('out',x.shape)
        return x  #50,32
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            prob = torch.sum(prob,dim=1).view(-1,1)
            return prob
        elif self.type =='cosine':
            prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            h = torch.cat([tf_embed, target_embed],dim=1)
            prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    def sample_indices(self,num_all_nodes, idx, num_samples=1000):
        """采样1000个节点，确保idx中256个节点包含在内"""
        # 将索引从tensor转换为numpy，如果它还不是
        idx = idx.cpu().numpy() if torch.is_tensor(idx) else idx
        
        idx = np.unique(idx)
      
        # 其余将要被采样的索引数
        num_other_samples = num_samples - len(idx)
        
        # 生成除了idx之外的所有节点索引
        other_indices = np.setdiff1d(np.arange(num_all_nodes), idx)

        # 从剩余节点中随机选取需要的数量来满足采样数
        other_sampled_indices = np.random.choice(other_indices, num_other_samples, replace=False)

        # 组合原始的idx和随机采样的结果
        sampled_indices = np.concatenate((idx, other_sampled_indices), axis=0)
        
        # 确保结果是打乱的
        np.random.shuffle(sampled_indices)

        return sampled_indices
    def resize(self, x_full, adj_full, sampled_idx):
        x_sampled = x_full[sampled_idx]
        #print(sampled_idx)
        adj_sampled = adj_full[sampled_idx][:, sampled_idx]
        return x_sampled, adj_sampled
    
# 假设有10000个节点，并且已经有了具体的idx

    def forward(self,x,adj,train_sample):
        #print(train_sample)
        adj = adj.to_dense()
        # 调整x和adj的大小以包含采样的数据
        sampled_indices = self.sample_indices(x.shape[0], train_sample)
        x_sampled, adj_sampled = self.resize(x, adj, sampled_indices)

        embed = self.encode(x_sampled,adj_sampled) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
               # 确保train_sample中的索引映射到采样后的节点上
        idx_map = {idx_value: i for i, idx_value in enumerate(sampled_indices)}
        #print(list(idx_map.keys()).count(7934))
        train_sampled_idx = torch.tensor([[idx_map[idx_i] for idx_i in idx_pair] for idx_pair in train_sample.tolist()], dtype=torch.long)
        
        # train_tf = tf_embed[train_sample[:,0]]
        # #ic(train_tf.shape)#train_tf 256 16
        # train_target = target_embed[train_sample[:, 1]]

        train_tf = tf_embed[train_sampled_idx[:, 0]]
        train_target = target_embed[train_sampled_idx[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output



class NeighborGCNGATLink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, num_head1, num_head2,
                alpha, device, type,reduction='mean'):
        super(NeighborGCNGATLink, self).__init__(device)
        self.device = device
        self.num_head1 = num_head1
        self.num_head2 = num_head2
        self.reduction = reduction
        self.type = type
        
        # 对于GCN，我们不使用多头机制
        self.gcn1 = GraphConvolution(input_dim, hidden1_dim)
        #self.gcn2 = GraphConvolution(hidden1_dim, hidden1_dim)

        # 对于GAT，我们根据头数创建多个注意力机制层
        # self.attention_heads = ModuleList([
        #     AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        # ])

        self.attention_heads = ModuleList([
            AttentionLayer(hidden1_dim, hidden2_dim, alpha) for _ in range(num_head2)
        ])
        # 如果采用拼接，更新hidden2_dim值为所有头输出拼接后的维度
        if reduction == 'concat':
            hidden2_dim = num_head2 * (hidden2_dim)
        elif reduction != 'mean':
            raise ValueError("Reduction method must be 'mean' or 'concat'.")

        self.hidden2_dim = hidden2_dim
        
        self.tf_linear1 = Linear(hidden2_dim, hidden3_dim)
        self.target_linear1 = Linear(hidden2_dim, hidden3_dim)

        self.tf_linear2 = Linear(hidden3_dim, output_dim)
        self.target_linear2 = Linear(hidden3_dim, output_dim)

        # 我们可以增加一个分类头，如果我们执行一个节点对分类任务
        self.linear = Linear(2 * output_dim, 2)
        
        self.reset_parameters()

    def reset_parameters(self):
        self.gcn1.reset_parameters()
        #self.gcn2.reset_parameters()
        for head in self.attention_heads:
            head.reset_parameters()

        # 初始化其余的全连接层
        nn.init.xavier_uniform_(self.tf_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear1.weight, gain=1.414)
        nn.init.xavier_uniform_(self.tf_linear2.weight, gain=1.414)
        nn.init.xavier_uniform_(self.target_linear2.weight, gain=1.414)
    

    def encode(self, x, adj):
        # GCN 编码

        x = F.relu(self.gcn1(x, adj))
        #x = F.relu(self.gcn2(x, adj))
        #ic('gcn',x.shape)
        # GAT 编码
        x_all = [att(x, adj) for att in self.attention_heads]
        #ic(x_all[1].shape)
        if self.reduction == 'concat':
            x = torch.cat(x_all, dim=1)
        elif self.reduction == 'mean':
            x = torch.mean(torch.stack(x_all, dim=0), dim=0)
        else:
            raise ValueError(f"Unknown reduction: {self.reduction}")
        #ic('out',x.shape)
        return x  #50,32
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            prob = torch.sum(prob,dim=1).view(-1,1)
            return prob
        elif self.type =='cosine':
            prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            h = torch.cat([tf_embed, target_embed],dim=1)
            prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def bfs_sample_neighbors(self,adj, gene, num_neighbors, include_self=True):
        import networkx as nx
        import scipy.sparse as sp
        sparse_tensor = adj.cpu()
        coo_indices = sparse_tensor._indices().numpy()
        coo_values = sparse_tensor._values().numpy()
        coo_shape = sparse_tensor.size()

        # 然后创建 scipy 的 COO 矩阵
        coo_matrix = sp.coo_matrix((coo_values, (coo_indices[0], coo_indices[1])), shape=coo_shape)

        # 下面你可以操作 scipy 的稀疏矩阵了，比如转换成 CSR 或 CSC 格式
        csr_matrix = coo_matrix.tocsr()
        nx_graph = nx.from_scipy_sparse_array(csr_matrix)
        
        gene = gene.cpu().numpy() if torch.is_tensor(gene) else gene
        genes = np.unique(gene).tolist()

        neighbors = nx.Graph()
        if include_self:
            for gene in genes:  
                neighbors.add_node(gene)
        bfs = nx.bfs_edges(nx_graph, gene)
        for u, v in bfs:
            if neighbors.number_of_nodes() == num_neighbors:
                break
            neighbors.add_node(v)

        for node in neighbors.nodes():
            for u, v, d in nx_graph.edges(node, data="weight"):
                if neighbors.has_node(u) and neighbors.has_node(v):
                    neighbors.add_weighted_edges_from([(u, v, d)])
        return neighbors
    
    def resize(self, x_full, adj_full, sampled_idx):
        x_sampled = x_full[sampled_idx]
        #print(sampled_idx)
        adj_sampled = adj_full[sampled_idx][:, sampled_idx]
        return x_sampled, adj_sampled
    
# 假设有10000个节点，并且已经有了具体的idx

    def forward(self,x,adj,train_sample):
        #print(train_sample)
        # 调整x和adj的大小以包含采样的数据
        #ic('input')
        neighbors = self.bfs_sample_neighbors(adj, train_sample,num_neighbors=2000)
        #ic('aaa')
        sampled_indices = list(neighbors.nodes)
        #ic('bbb')
        adj = adj.to_dense()
        #ic('ccc')
        x_sampled, adj_sampled = self.resize(x, adj, sampled_indices)

        #ic('ddd')
        embed = self.encode(x_sampled,adj_sampled) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
            # 确保train_sample中的索引映射到采样后的节点上
        idx_map = {idx_value: i for i, idx_value in enumerate(sampled_indices)}
        #print(list(idx_map.keys()).count(7934))
        train_sampled_idx = torch.tensor([[idx_map[idx_i] for idx_i in idx_pair] for idx_pair in train_sample.tolist()], dtype=torch.long)
        
        # train_tf = tf_embed[train_sample[:,0]]
        # #ic(train_tf.shape)#train_tf 256 16
        # train_target = target_embed[train_sample[:, 1]]

        train_tf = tf_embed[train_sampled_idx[:, 0]]
        train_target = target_embed[train_sampled_idx[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output

class GraphSAGELink(BaseModel):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training, num_head2,
                 flag, device, type,hidden_dim,visualize,reduction='mean'):#
        super(GraphSAGELink, self).__init__(device)
        from torch_geometric.nn import SAGEConv
        self.device = device
        #self.num_head1 = num_head1
        self.num_head2 = num_head2
        self.reduction = reduction
        self.type = type
        self.flag = flag
        
        # GraphSAGE层
        self.sage1 = SAGEConv(input_dim, hidden1_dim)
        self.sage2 = SAGEConv(hidden1_dim, hidden2_dim)

        self.hidden2_dim = hidden2_dim
        
        self.tf_linear1 = Linear(hidden2_dim, hidden3_dim)
        self.target_linear1 = Linear(hidden2_dim, hidden3_dim)

        self.tf_linear2 = Linear(hidden3_dim, output_dim)
        self.target_linear2 = Linear(hidden3_dim, output_dim)

        # 我们可以增加一个分类头，如果我们执行一个节点对分类任务
        self.linear = Linear(2 * output_dim, 2)
        self.decode_linear = Linear(output_dim, 2)
        self.reset_parameters()

    def reset_parameters(self):
        self.sage1.reset_parameters()
        self.sage2.reset_parameters()
        self.tf_linear1.reset_parameters()
        self.target_linear1.reset_parameters()
        self.tf_linear2.reset_parameters()
        self.target_linear2.reset_parameters()
        
    def encode(self, x, adj):
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()

            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()

            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)

        x = F.relu(self.sage1(x, edge_index))
        x = F.dropout(x, p=0.01, training=self.training)
        x = F.relu(self.sage2(x, edge_index))
        return x  #50,32
    
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
        
    def forward(self,x,adj,train_sample,gene_num):
        #ic(x.shape,adj.shape,train_sample.shape)#x 910,758  adj 910,910  train_sample 256 2
        #ic('-----',x.shape)
        embed = self.encode(x,adj) 
        #ic(embed.shape) #embed 910 64
        #ic(embed.shape)
        tf_embed = self.tf_linear1(embed)
        #ic(tf_embed.shape)#tf_embed 910 32
        tf_embed = F.leaky_relu(tf_embed)
        tf_embed = F.dropout(tf_embed,p=0.01)
        tf_embed = self.tf_linear2(tf_embed)#tf_embed 910 16
        #ic(tf_embed.shape)
        tf_embed = F.leaky_relu(tf_embed)
        target_embed = self.target_linear1(embed)
        #ic(target_embed.shape)#target_embed 910 32
        target_embed = F.leaky_relu(target_embed)
        target_embed = F.dropout(target_embed, p=0.01)
        target_embed = self.target_linear2(target_embed)#target_embed 910 16
        #ic(target_embed.shape)
        target_embed = F.leaky_relu(target_embed)
        self.tf_ouput = tf_embed
        self.target_output = target_embed
        train_tf = tf_embed[train_sample[:,0]]
        #ic(train_tf.shape)#train_tf 256 16
        train_target = target_embed[train_sample[:, 1]]
        #ic(train_target.shape)#train_target 256 16
        pred = self.decode(train_tf, train_target)
        #ic(pred.shape)#pred 256 1
        return pred

    def get_embedding(self):
        return self.tf_ouput, self.target_output

from torch.nn import Sequential, Linear
from torch_geometric.nn import SAGEConv, GATConv
from torch_geometric.utils import to_dense_adj, dense_to_sparse

class PinSageLink(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training, num_head2,
                 flag, device, type,reduction='mean'):
        super(PinSageLink, self).__init__()
        self.device = device
        self.type = type
        # 定义两个PinSage卷积层
        self.conv1 = SAGEConv(input_dim, hidden1_dim, normalize=True, aggr='mean')
        self.conv2 = SAGEConv(hidden1_dim, output_dim, normalize=True, aggr='mean')
        self.flag = flag
        # 在类初始化中定义采样和池化策略

        #self.num_neighbors = num_neighbors

    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj, train_sample):
        #-----------跑起来慢------------------------------
        # if adj.is_sparse:
        #     adj = adj.to_dense()
        # if adj.is_cuda:
        #     adj = adj.cpu()
        # edge_index = adj.nonzero().t().contiguous()
        # edge_index = edge_index.to(self.device)  # 跑起来很慢

        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()

            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()

            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
            
        # 第一层PinSage卷积
        x = self.conv1(x, edge_index)
        # 应用非线性激活函数和dropout
        x = F.relu(x)
        x = F.dropout(x, p=0.5, training=self.training)

        # 第二层PinSage卷积
        x = self.conv2(x, edge_index)
        x = F.relu(x)

        # 分别获取训练节点对中的两部分节点的嵌入
        tf_embed = x[train_sample[:,0]]
        target_embed = x[train_sample[:,1]]

        # 解码器可以保留前面的
        pred = self.decode(tf_embed, target_embed)
        return pred

from torch_geometric.utils import add_self_loops, degree
from torch_geometric.nn.inits import glorot, zeros
from torch_scatter import scatter_add

class FastGCNConv(nn.Module):
    def __init__(self, in_channels, out_channels, bias=True):
        super(FastGCNConv, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        
        self.weight = nn.Parameter(torch.Tensor(in_channels, out_channels))
        if bias:
            self.bias = nn.Parameter(torch.Tensor(out_channels))
        else:
            self.register_parameter('bias', None)
        self.reset_parameters()

    def reset_parameters(self):
        glorot(self.weight)
        zeros(self.bias)

    def forward(self, x, edge_index, edge_weight=None):
        x = torch.matmul(x, self.weight)

        if edge_weight is None:
            edge_index, _ = add_self_loops(edge_index, num_nodes=x.size(0))
            row, col = edge_index
            deg = scatter_add(torch.ones_like(col), row, dim=0, dim_size=x.size(0))
        else:
            edge_index, edge_weight = add_self_loops(edge_index, edge_weight, fill_value=1, num_nodes=x.size(0))
            row, col = edge_index
            deg = scatter_add(edge_weight, row, dim=0, dim_size=x.size(0))

        deg_inv_sqrt = deg.pow(-0.5)
        deg_inv_sqrt[deg_inv_sqrt == float('inf')] = 0

        if edge_weight is None:
            edge_weight = torch.ones_like(row, dtype=x.dtype)

        edge_weight = deg_inv_sqrt[row] * edge_weight * deg_inv_sqrt[col]
        out = scatter_add(x[row] * edge_weight.view(-1, 1), col, dim=0, dim_size=x.size(0))

        if self.bias is not None:
            out += self.bias

        return out

class FastGCNLink(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training, num_head2,
                 flag, device, type,reduction='mean'):
        super(FastGCNLink, self).__init__()
        # 使用 FastGCN 卷积层代替 PinSage 层
        self.flag = flag
        self.device = device
        self.conv1 = FastGCNConv(input_dim, hidden1_dim)
        self.conv2 = FastGCNConv(hidden1_dim, hidden2_dim)
        self.conv3 = FastGCNConv(hidden2_dim, output_dim)
        self.decode_linear = Linear(output_dim, 2)
        self.type = type
        if self.type == 'MLP':
            # MLP 解码器可能包含多个层
            self.linear = nn.Sequential(
                nn.Linear(2*output_dim, hidden2_dim),
                nn.ReLU(),
                nn.Linear(hidden2_dim, 1)
            )

    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
        
    def forward(self, x, adj, train_sample):
        # FastGCN也需要在整个图上运行卷积操作，因此需要传入完整的edge_index，但是对应权重的计算可能需要适配
        # 书写卷积层的前向传播
         #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()

            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()

            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)

        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.conv2(x, edge_index))
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv3(x, edge_index)

        # 获取train_sample中对应的节点嵌入
        tf_embed = x[train_sample[:, 0], :]
        target_embed = x[train_sample[:, 1], :]

        # 调用解码器（具体取决于自己的解码器实现）
        pred = self.decode(tf_embed, target_embed)
        return pred
from torch_geometric.nn import GCNConv, ClusterGCNConv

class ClusterGCNLink(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training, num_head2,
                 flag, device, type,reduction='mean'):
        super(ClusterGCNLink, self).__init__()
        # 不需要 flag，因为 ClusterGCN 假设是对无向图进行操作
        self.device = device
        self.conv1 = ClusterGCNConv(input_dim, hidden1_dim)
        self.conv2 = ClusterGCNConv(hidden1_dim, hidden2_dim)
        self.conv3 = ClusterGCNConv(hidden2_dim, output_dim)
        self.decode_linear = Linear(output_dim, 2)
        self.flag = flag
        self.type = type
        if self.type == 'MLP':
            self.linear = torch.nn.Sequential(
                torch.nn.Linear(2 * output_dim, hidden2_dim),
                torch.nn.ReLU(),
                torch.nn.Linear(hidden2_dim, 1)
            )
        # 省略了原模型中解码器部分的定义，因为它和模型类型的转换无关
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
        
    def forward(self, x,adj, train_sample):

        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()

            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()

            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        # ClusterGCN 不需要单独的邻接矩阵，因为它希望在每次前向传播时对子图进行操作
        # 运行 ClusterGCN 卷积操作
        x = F.relu(self.conv1(x, edge_index))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.conv2(x, edge_index))
        x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv3(x, edge_index)

        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x[train_sample[:, 0], :]
        target_embed = x[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        return pred

#from torch_geometric.nn import TransformerConv
from models.weight_transformer_conv import TransformerConv
class GraphTransformerLink(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type,hidden_dim,visualize,reduction):
        super(GraphTransformerLink, self).__init__()
        dropout=0.1
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize

        # 图 Transformer 卷积层
        self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=5, dropout=dropout, edge_dim=edge_dim)
        self.conv2 = TransformerConv(hidden1_dim * 5, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        #self.conv3 = TransformerConv(hidden2_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv3 = TransformerConv(hidden2_dim*3*1, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,gene_num):
        edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        
        x1 = F.relu(self.conv1(x, edge_index, edge_attr))
        x = F.dropout(x1, p=0.5, training=self.training)
        x2 = F.relu(self.conv2(x, edge_index, edge_attr))
        x = F.dropout(x2, p=0.5, training=self.training)
        x3 = self.conv3(x, edge_index, edge_attr)
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize:
            return pred,tf_embed,target_embed,x1,self.conv1.alpha,x2,self.conv2.alpha,x3,self.conv3.alpha
        else:
            return pred

class GraphNodeFeature(nn.Module):
    """
    Compute node features for each node in the graph.
    """

    def __init__(
        self,hidden_dim=1134, num_in_degree=86+1, num_out_degree=7485+1,
    ): #japonica hidden_dim=1134,num_in_degree=81+1, num_out_degree=6203+1,
        #indica  num_in_degree=86+1, num_out_degree=7485+1,
        super(GraphNodeFeature, self).__init__()
        self.hidden_dim =hidden_dim
        self.in_degree_encoder = nn.Embedding(num_embeddings=num_in_degree, embedding_dim=hidden_dim, padding_idx=0)
        self.out_degree_encoder = nn.Embedding(num_embeddings=num_out_degree, embedding_dim=hidden_dim, padding_idx=0)
        
        self.in_degree = None
        self.out_degree = None
        
    def forward(self,x,adj):
        #print(x.max(),x.min())
        if self.in_degree is None:
            adj_dense = adj.to_dense()
            self.out_degree = torch.sum(adj_dense, dim=1).long()
            #max_out_degree = torch.max(self.out_degree)
            self.in_degree = torch.sum(adj_dense, dim=0).long()
            #max_in_degree = torch.max(self.in_degree)
        #self.hidden_dim = hidden_dim
        #max_in_degree = torch.max(in_degree)
        #ic(max_out_degree,max_in_degree)
        #ic(x.shape,self.in_degree.shape,self.out_degree.shape)
        node_feature = x + self.in_degree_encoder(self.in_degree) + self.out_degree_encoder(self.out_degree)
        #print(x.dtype,node_feature.dtype,x.shape,node_feature.shape)
       
        #print(node_feature.max(),node_feature.min())
        return node_feature

class GraphAttnBias(nn.Module):
    """
    Compute attention bias for each head.
    """

    def __init__(self, num_heads, num_spatial):
        super(GraphAttnBias, self).__init__()
        self.num_heads = num_heads
        self.spatial_pos_encoder = nn.Embedding(num_spatial, num_heads, padding_idx=0)
        self.graph_token_virtual_distance = nn.Embedding(1, num_heads)

    def forward(self, x, spatial_pos, attn_bias):
        # 确保 spatial_pos 是正确的类型
        spatial_pos = spatial_pos.long()  # 转换索引为长整型
        spatial_pos_bias = self.spatial_pos_encoder(spatial_pos).permute(0, 3, 1, 2)
        attn_bias = attn_bias.unsqueeze(1)  # 最小化扩展操作
        graph_attn_bias = attn_bias + spatial_pos_bias  # 减少重复和过多的扩展
        return graph_attn_bias

from models.weight_transformer_conv import TransformerConv    
class Graphormer1Link(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, hidden_dim,visualize,reduction,**kwgs):
        super(Graphormer1Link, self).__init__()
        dropout=0.1
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize
        adj = kwgs.get('adj',None)
        adj_dense = adj.to_dense()
        out_degree = torch.sum(adj_dense, dim=1).long()
        max_out_degree = torch.max(out_degree)
        in_degree = torch.sum(adj_dense, dim=0).long()
        max_in_degree = torch.max(in_degree)    
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=hidden_dim,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        # 图 Transformer 卷积层

        # self.conv1 = SAGEConv(input_dim, hidden1_dim)
        # self.conv2 = SAGEConv(hidden1_dim, hidden2_dim)
        self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        # self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,gene_num):
        edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            if not adj.is_sparse:
                adj = adj.to_sparse()
            adj_matrix_coalesced = adj.coalesce()
            #adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(x.shape)
        x = self.GraphNodeFeature(x,adj)
        x0=x
        #ic(x.shape,edge_index.shape)
        x1 = F.relu(self.conv1(x, edge_index, edge_attr))
        x = F.dropout(x1, p=0.5, training=self.training)
        x2 = F.relu(self.conv2(x, edge_index, edge_attr))
        x = F.dropout(x2, p=0.5, training=self.training)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)
        x3 = self.conv3(x, edge_index, edge_attr)
        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize: #self.conv1.alpha,self.conv2.alpha,,self.conv3.alpha
            return x0,edge_index,pred,tf_embed,target_embed,x1,x2,x3
        else:
            return pred

from models.weight_transformer_conv import TransformerConv    
class Graphormer1Link_type(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, hidden_dim,visualize,reduction,**kwgs):
        super(Graphormer1Link_type, self).__init__()
        dropout=0.1
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize
        adj = kwgs.get('adj',None)
        adj_dense = adj.to_dense()
        out_degree = torch.sum(adj_dense, dim=1).long()
        max_out_degree = torch.max(out_degree)
        in_degree = torch.sum(adj_dense, dim=0).long()
        max_in_degree = torch.max(in_degree)    
        self.layer_options = {
            'japonica_CK': nn.Linear(604, 1000).to(device),
            'japonica_CD': nn.Linear(636, 1000).to(device),
            'japonica_CE': nn.Linear(85, 1000).to(device),
            'indica_CK': nn.Linear(333, 1000).to(device),
            'indica_CD': nn.Linear(265, 1000).to(device),
            'indica_CE': nn.Linear(56, 1000).to(device),
            'zeamays_CK': nn.Linear(745, 1000).to(device),
            'zeamays_CD': nn.Linear(780, 1000).to(device),
            'zeamays_CE': nn.Linear(10, 1000).to(device),
        }
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=1000,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        # 图 Transformer 卷积层

        # self.conv1 = SAGEConv(input_dim, hidden1_dim)
        # self.conv2 = SAGEConv(hidden1_dim, hidden2_dim)#input_dim
        self.conv1 = TransformerConv(1000, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        # self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,gene_num):

        if x.shape[1] == 604:
            x = self.layer_options['japonica_CK'](x)
        elif x.shape[1] == 636:
            x = self.layer_options['japonica_CD'](x)
        elif x.shape[1] == 85:
            x = self.layer_options['japonica_CE'](x)
        elif x.shape[1] == 333:
            x = self.layer_options['indica_CK'](x)
        elif x.shape[1] == 265:
            x = self.layer_options['indica_CD'](x)
        elif x.shape[1] == 56:
            x = self.layer_options['indica_CE'](x)
        elif x.shape[1] == 745:
            x = self.layer_options['zeamays_CK'](x)
        elif x.shape[1] == 780:
           # print(self.layer_options['zeamays_CD'])
            x = self.layer_options['zeamays_CD'](x)
        elif x.shape[1] == 10:
            x = self.layer_options['zeamays_CE'](x)

        edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            if not adj.is_sparse:
                adj = adj.to_sparse()
            adj_matrix_coalesced = adj.coalesce()
            #adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(x.shape)
        x = self.GraphNodeFeature(x,adj)
        x0=x
        #ic(x.shape,edge_index.shape)
        x1 = F.relu(self.conv1(x, edge_index, edge_attr))
        x = F.dropout(x1, p=0.5, training=self.training)
        x2 = F.relu(self.conv2(x, edge_index, edge_attr))
        x = F.dropout(x2, p=0.5, training=self.training)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)
        x3 = self.conv3(x, edge_index, edge_attr)
        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize: #self.conv1.alpha,self.conv2.alpha,,self.conv3.alpha
            return adj,x0,edge_index,pred,tf_embed,target_embed,x1,x2,x3
        else:
            return pred

class Graphormer1Link_gate(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, hidden_dim,visualize,reduction,**kwgs):
        super(Graphormer1Link_gate, self).__init__()
        dropout=0.1
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize
        adj = kwgs.get('adj',None)
        adj_dense = adj.to_dense()
        out_degree = torch.sum(adj_dense, dim=1).long()
        max_out_degree = torch.max(out_degree)
        in_degree = torch.sum(adj_dense, dim=0).long()
        max_in_degree = torch.max(in_degree)    
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=hidden_dim,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        # 图 Transformer 卷积层

        # self.conv1 = SAGEConv(input_dim, hidden1_dim)
        # self.conv2 = SAGEConv(hidden1_dim, hidden2_dim)
        self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.norm1 = nn.LayerNorm(hidden1_dim*3)
        self.gate1 = nn.Linear(input_dim, hidden1_dim*3)
        self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.norm2 = nn.LayerNorm(hidden2_dim*3)
        self.gate2 = nn.Linear(hidden1_dim * 3, hidden2_dim*3)
        self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)
        self.norm3 = nn.LayerNorm(output_dim*1)
        self.gate3 = nn.Linear(hidden2_dim*3, output_dim)

        # self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,gene_num):
        edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            if not adj.is_sparse:
                adj = adj.to_sparse()
            adj_matrix_coalesced = adj.coalesce()
            #adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(x.shape)
        x = self.GraphNodeFeature(x,adj)
        x0=x
        #ic(x.shape,edge_index.shape)
        x1 = F.relu(self.norm1(self.conv1(x, edge_index, edge_attr)))
        #ic(x1.shape,self.gate1(x0).shape)
        x1=x1+self.gate1(x0)*F.sigmoid(x1)
        x = F.dropout(x1, p=0.5, training=self.training)
        x2 = F.relu(self.norm2(self.conv2(x, edge_index, edge_attr)))
        x2=x2+self.gate2(x)*F.sigmoid(x2)
        x = F.dropout(x2, p=0.5, training=self.training)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)
        x3 = F.relu(self.norm3(self.conv3(x, edge_index, edge_attr)))
        x3=x3+self.gate3(x)*F.sigmoid(x3)
        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize: #self.conv1.alpha,self.conv2.alpha,,self.conv3.alpha
            return adj,x0,edge_index,pred,tf_embed,target_embed,x1,x2,x3
        else:
            return x0, pred



from models.GraphormerconvLayer import GraphormerConv
class Graphormer2Link_gate(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim,
                flag, device, type, hidden_dim,visualize,**kwgs):
        super(Graphormer2Link_gate, self).__init__()
        dropout=0.1
        #edge_dim=hidden1_dim
        self.flag = flag
        self.device = device
        #self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.hidden1_dim=hidden1_dim
        self.hidden2_dim=hidden2_dim
        self.output_dim=output_dim

        self.visualize = visualize
        # adj = kwgs.get('adj',None)
        # adj_dense = adj.to_dense()
        # out_degree = torch.sum(adj_dense, dim=1).long()
        # max_out_degree = torch.max(out_degree)
        # in_degree = torch.sum(adj_dense, dim=0).long()
        # max_in_degree = torch.max(in_degree)    
        max_out_degree = kwgs.get('max_out_degree',None)
        max_in_degree = kwgs.get('max_in_degree',None)

        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=hidden_dim,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        self.GraphAttnBias = GraphAttnBias(num_heads=3,num_spatial=7+1)
        # 图 Transformer 卷积层
        self.attn_bias1=nn.Linear(1,3)
        self.transform1 = nn.Linear(input_dim, hidden1_dim)
        self.conv1 = GraphormerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=hidden1_dim)
        self.norm1 = nn.LayerNorm(hidden1_dim*3)
        self.gate1 = nn.Linear(input_dim, hidden1_dim*3)

        self.attn_bias2=nn.Linear(1,3)
        self.transform2 = nn.Linear(hidden1_dim*3 , hidden2_dim)
        self.conv2 = GraphormerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=hidden2_dim)
        self.norm2 = nn.LayerNorm(hidden2_dim*3)
        self.gate2 = nn.Linear(hidden1_dim * 3, hidden2_dim*3)
        #self.conv3 = TransformerConv(hidden2_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)

        self.attn_bias3=nn.Linear(1,1)
        self.transform3 = nn.Linear(hidden2_dim*3, output_dim)
        self.conv3 = GraphormerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=output_dim,concat=False)
        self.norm3 = nn.LayerNorm(output_dim*1)
        self.gate3 = nn.Linear(hidden2_dim*3, output_dim)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,att_bias):

        #ic(space_encoding_tensor.max())
        
       
        # if self.flag:
        #     prob = self.decode_linear(prob)
        # else:
        #     prob = torch.sum(prob,dim=1).view(-1,1)
        # edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(space_encoding_tensor)

        rows = edge_index[0]  # 从 edge_index 第一行获取起始节点索引
        cols = edge_index[1]  # 从 edge_index 第二行获取结束节点索引

        # 使用这些行和列索引从 space_encoding_tensor 中提取对应的值
        selected_values = att_bias[rows, cols]

        # 如果您需要保持 selected_values 的形状为 [num_edges, 1] 或类似格式，可以调整形状
        att_bias = selected_values#.unsqueeze(1)
        att_bias = att_bias.unsqueeze(1)
        #ic(att_bias.shape)
        #att_bias = None#self.GraphAttnBias(x,space_encoding_tensor,att_bias)
        #ic(x.shape,adj.shape)
        x = self.GraphNodeFeature(x,adj)
        #ic(x.shape, edge_index.shape, edge_attr.shape)
        x0 = x

        att_bias1 = self.attn_bias1(att_bias)
        edge_attr1 = self.transform1(x[edge_index[0]] + x[edge_index[1]])
        x1 = F.relu(self.norm1(self.conv1(x, edge_index, edge_attr1,att_bias=att_bias1)))
        x1=x1+self.gate1(x0)*F.sigmoid(x1)
        x = F.dropout(x1, p=0.5)


        att_bias2 = self.attn_bias2(att_bias)
        edge_attr2 = self.transform2(x[edge_index[0]] + x[edge_index[1]])
        x2 = F.relu(self.norm2(self.conv2(x, edge_index, edge_attr2,att_bias=att_bias2)))
        x2=x2+self.gate2(x)*F.sigmoid(x2)
        x = F.dropout(x2, p=0.5)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)


        att_bias3 = self.attn_bias3(att_bias)
        edge_attr3 = self.transform3(x[edge_index[0]] + x[edge_index[1]])
        x3 = F.relu(self.norm3(self.conv3(x, edge_index, edge_attr3,att_bias=att_bias3)))
        x3=x3+self.gate3(x)*F.sigmoid(x3)
        # 从节点特征表示中抽取特定训练样本的特征
        #ic(x3.shape)
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]
        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)

        #ic(pred.shape)
        if self.visualize:
            return adj,x0,edge_index,pred,tf_embed,target_embed,x1,self.conv1.alpha,x2,self.conv2.alpha,x3,self.conv3.alpha
        else:
            return pred


from models.GraphormerconvLayer import GraphormerConv
class Graphormer2Link_gate_type(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim,
                flag, device, type, hidden_dim,visualize,**kwgs):
        super(Graphormer2Link_gate_type, self).__init__()
        dropout=0.1
        #edge_dim=hidden1_dim
        self.flag = flag
        self.device = device
        #self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.hidden1_dim=hidden1_dim
        self.hidden2_dim=hidden2_dim
        self.output_dim=output_dim

        self.visualize = visualize
        #
        self.layer_options =nn.ModuleDict({
            'japonica_CK': nn.Linear(602, 1000).to(device),
            'japonica_CD': nn.Linear(635, 1000).to(device),
            'japonica_CE': nn.Linear(85, 1000).to(device),
            'japonica_all': nn.Linear(1415, 1000).to(device),
            'indica_CK': nn.Linear(333, 1000).to(device),
            'indica_CD': nn.Linear(263, 1000).to(device),
            'indica_CE': nn.Linear(56, 1000).to(device),
            'zeamays_CK': nn.Linear(743, 1000).to(device),
            'zeamays_CD': nn.Linear(777, 1000).to(device),
            'zeamays_CE': nn.Linear(10, 1000).to(device),
            # 'wheat_CK': nn.Linear(256, 1000).to(device),
            # 'wheat_CD': nn.Linear(452, 1000).to(device),
            # 'wheat_CE': nn.Linear(21, 1000).to(device),
            # 'sorghum_CK': nn.Linear(194, 1000).to(device),
            # 'sorghum_CD': nn.Linear(265, 1000).to(device),
            # 'sorghum_CE': nn.Linear(5, 1000).to(device),
            # 'homo_C3_CK': nn.Linear(327, 1000).to(device),
            # 'homo_C3_CD': nn.Linear(257, 1000).to(device),
            # 'homo_C3_CE': nn.Linear(56, 1000).to(device),
            # 'homo_C4_CK': nn.Linear(937, 1000).to(device),
            # 'homo_C4_CD': nn.Linear(1042, 1000).to(device),
            # 'homo_C4_CE': nn.Linear(10, 1000).to(device),
            # 'homo_CK': nn.Linear(937, 1000).to(device),
            # 'homo_CD': nn.Linear(1042, 1000).to(device),
            # 'homo_CE': nn.Linear(10, 1000).to(device),
            'homo_CK': nn.Linear(1676, 1000).to(device),
            'homo_CD': nn.Linear(1675, 1000).to(device),
            'homo_CE': nn.Linear(151, 1000).to(device),
        })
        #self.layer = None
        self.max_out_degree = kwgs.get('max_out_degree',None)
        self.max_in_degree = kwgs.get('max_in_degree',None)
        #ic(self.max_out_degree,self.max_in_degree)
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=1000,num_in_degree=self.max_in_degree+1,num_out_degree=self.max_out_degree+1)
        
        self.GraphAttnBias = GraphAttnBias(num_heads=3,num_spatial=7+1)
        # 图 Transformer 卷积层
        self.attn_bias1=nn.Linear(1,3)
        self.transform1 = nn.Linear(1000, hidden1_dim)
        self.conv1 = GraphormerConv(1000, hidden1_dim, heads=3, dropout=dropout, edge_dim=hidden1_dim)
        self.norm1 = nn.LayerNorm(hidden1_dim*3)
        self.gate1 = nn.Linear(1000, hidden1_dim*3)

        self.attn_bias2=nn.Linear(1,3)
        self.transform2 = nn.Linear(hidden1_dim*3 , hidden2_dim)
        self.conv2 = GraphormerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=hidden2_dim)
        self.norm2 = nn.LayerNorm(hidden2_dim*3)
        self.gate2 = nn.Linear(hidden1_dim * 3, hidden2_dim*3)
        #self.conv3 = TransformerConv(hidden2_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)

        self.attn_bias3=nn.Linear(1,1)
        self.transform3 = nn.Linear(hidden2_dim*3, output_dim)
        self.conv3 = GraphormerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=output_dim,concat=False)
        self.norm3 = nn.LayerNorm(output_dim*1)
        self.gate3 = nn.Linear(hidden2_dim*3, output_dim)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
    
    # def init_adj(self,adj):
    #     #adj = kwgs.get('adj',None)
    #     adj_dense = adj.to_dense()
    #     out_degree = torch.sum(adj_dense, dim=1).long()
    #     self.max_out_degree = torch.max(out_degree)
    #     in_degree = torch.sum(adj_dense, dim=0).long()
    #     self.max_in_degree = torch.max(in_degree)    


    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,att_bias):
        #print(x.shape)
        input=x
        #print(x.device,self.layer_options['japonica_CK'].weight.device)
        if x.shape[1] == 602:
            x = self.layer_options['japonica_CK'](x)
        elif x.shape[1] == 635:
            x = self.layer_options['japonica_CD'](x)
        elif x.shape[1] == 85:
            x = self.layer_options['japonica_CE'](x)
        elif x.shape[1] == 1415:
            x = self.layer_options['japonica_all'](x)
        elif x.shape[1] == 333:
            x = self.layer_options['indica_CK'](x)
        elif x.shape[1] == 263:
            x = self.layer_options['indica_CD'](x)
        elif x.shape[1] == 56:
            x = self.layer_options['indica_CE'](x)
        elif x.shape[1] == 743:
            x = self.layer_options['zeamays_CK'](x)
        elif x.shape[1] == 777:
            x = self.layer_options['zeamays_CD'](x)
        elif x.shape[1] == 10:
            x = self.layer_options['zeamays_CE'](x)
        # elif x.shape[1] == 256:
        #     x = self.layer_options['wheat_CK'](x)
        # elif x.shape[1] == 452:
        #     x = self.layer_options['wheat_CD'](x)
        # elif x.shape[1] == 21:
        #     x = self.layer_options['wheat_CE'](x)
        # elif x.shape[1] == 194:
        #     x = self.layer_options['sorghum_CK'](x)
        # elif x.shape[1] == 265:
        #     x = self.layer_options['sorghum_CD'](x)
        # elif x.shape[1] == 5:
        #     x = self.layer_options['sorghum_CE'](x)
        # elif x.shape[1] == 327:
        #     x = self.layer_options['homo_C3_CK'](x)
        # elif x.shape[1] == 257:
        #     x = self.layer_options['homo_C3_CD'](x)
        # elif x.shape[1] == 56:
        #     x = self.layer_options['homo_C3_CE'](x)
        # elif x.shape[1] == 937:
        #     x = self.layer_options['homo_C4_CK'](x)
        # elif x.shape[1] == 1042:
        #    # print(self.layer_options['zeamays_CD'])
        #     x = self.layer_options['homo_C4_CD'](x)
        # elif x.shape[1] == 10:
        #     x = self.layer_options['homo_C4_CE'](x)
        elif x.shape[1] == 1676:#937:
            x = self.layer_options['homo_CK'](x)
        elif x.shape[1] == 1675:#1042:
           # print(self.layer_options['zeamays_CD'])
            x = self.layer_options['homo_CD'](x)
        elif x.shape[1] == 151:#10:
            x = self.layer_options['homo_CE'](x)
        #ic(x)
        #print(x.shape)
        # input_dim = x.shape[1]
        # self.layer = nn.Linear(input_dim, 1000).to(x.device)
        # x=self.layer(x)

        #ic(space_encoding_tensor.max())
        #print(adj)
       
        # if self.flag:
        #     prob = self.decode_linear(prob)
        # else:
        #     prob = torch.sum(prob,dim=1).view(-1,1)
        # edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(space_encoding_tensor)
        #print(x.shape,edge_index.shape,edge_index.max())
        rows = edge_index[0]  # 从 edge_index 第一行获取起始节点索引
        cols = edge_index[1]  # 从 edge_index 第二行获取结束节点索引
        #ic(rows.max(),cols.max())
        # 使用这些行和列索引从 space_encoding_tensor 中提取对应的值
        selected_values = att_bias[rows, cols]

        # 如果您需要保持 selected_values 的形状为 [num_edges, 1] 或类似格式，可以调整形状
        att_bias = selected_values#.unsqueeze(1)
        att_bias = att_bias.unsqueeze(1)
        #ic(att_bias.shape)
        #att_bias = None#self.GraphAttnBias(x,space_encoding_tensor,att_bias)
        #ic(x)
        x = self.GraphNodeFeature(x,adj)
        #ic(x)
        #ic(x.shape, edge_index.shape, edge_attr.shape)
        x0 = x

        att_bias1 = self.attn_bias1(att_bias)
        #ic(x.shape,edge_index.shape)
        edge_attr1 = self.transform1(x[edge_index[0]] + x[edge_index[1]])
        x1 = F.relu(self.norm1(self.conv1(x, edge_index, edge_attr1,att_bias=att_bias1)))
        x1=x1+self.gate1(x0)*F.sigmoid(x1)
        x = F.dropout(x1, p=0.5)

        att_bias2 = self.attn_bias2(att_bias)
        edge_attr2 = self.transform2(x[edge_index[0]] + x[edge_index[1]])
        x2 = F.relu(self.norm2(self.conv2(x, edge_index, edge_attr2,att_bias=att_bias2)))
        x2=x2+self.gate2(x)*F.sigmoid(x2)
        x = F.dropout(x2, p=0.5)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)

        att_bias3 = self.attn_bias3(att_bias)
        edge_attr3 = self.transform3(x[edge_index[0]] + x[edge_index[1]])
        x3 = F.relu(self.norm3(self.conv3(x, edge_index, edge_attr3,att_bias=att_bias3)))
        x3=x3+self.gate3(x)*F.sigmoid(x3)
        # 从节点特征表示中抽取特定训练样本的特征
        #ic(x3,train_sample)
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]
        
        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)

        
        if self.visualize:
            return input,adj,x0,edge_index,pred,tf_embed,target_embed,x1,self.conv1.alpha,x2,self.conv2.alpha,x3,self.conv3.alpha
        else:
            return pred

class Graphormer2Link_gate_type_small(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim,
                flag, device, type, hidden_dim,visualize,**kwgs):
        super(Graphormer2Link_gate_type_small, self).__init__()
        dropout=0.1
        #edge_dim=hidden1_dim
        self.flag = flag
        self.device = device
        #self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.hidden1_dim=hidden1_dim
        self.hidden2_dim=hidden2_dim
        self.output_dim=output_dim

        self.visualize = visualize

        self.layer_options =nn.ModuleDict({
            'japonica_CK': nn.Linear(602, 1000).to(device),
            'japonica_CD': nn.Linear(635, 1000).to(device),
            'japonica_CE': nn.Linear(85, 1000).to(device),
            'indica_CK': nn.Linear(333, 1000).to(device),
            'indica_CD': nn.Linear(263, 1000).to(device),
            'indica_CE': nn.Linear(56, 1000).to(device),
            'zeamays_CK': nn.Linear(743, 1000).to(device),
            'zeamays_CD': nn.Linear(777, 1000).to(device),
            'zeamays_CE': nn.Linear(10, 1000).to(device),
            'homo_CK': nn.Linear(1676, 1000).to(device),
            'homo_CD': nn.Linear(1675, 1000).to(device),
            'homo_CE': nn.Linear(151, 1000).to(device),
        })
        #self.layer = None
        self.max_out_degree = kwgs.get('max_out_degree',None)
        self.max_in_degree = kwgs.get('max_in_degree',None)
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=1000,num_in_degree=self.max_in_degree+1,num_out_degree=self.max_out_degree+1)
        
        self.GraphAttnBias = GraphAttnBias(num_heads=3,num_spatial=7+1)
        # 图 Transformer 卷积层
        self.attn_bias1=nn.Linear(1,3)
        self.transform1 = nn.Linear(1000, hidden1_dim)
        self.conv1 = GraphormerConv(1000, hidden1_dim, heads=3, dropout=dropout, edge_dim=hidden1_dim)
        self.norm1 = nn.LayerNorm(hidden1_dim*3)
        self.gate1 = nn.Linear(1000, hidden1_dim*3)

        # self.attn_bias2=nn.Linear(1,3)
        # self.transform2 = nn.Linear(hidden1_dim*3 , hidden2_dim)
        # self.conv2 = GraphormerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=hidden2_dim)
        # self.norm2 = nn.LayerNorm(hidden2_dim*3)
        # self.gate2 = nn.Linear(hidden1_dim * 3, hidden2_dim*3)
        #self.conv3 = TransformerConv(hidden2_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)

        self.attn_bias3=nn.Linear(1,1)
        self.transform3 = nn.Linear(hidden2_dim*3, output_dim)
        self.conv3 = GraphormerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=output_dim,concat=False)
        self.norm3 = nn.LayerNorm(output_dim*1)
        self.gate3 = nn.Linear(hidden2_dim*3, output_dim)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )

    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,att_bias):
        #print(x.shape)
        input=x
        #print(x.device,self.layer_options['japonica_CK'].weight.device)
        if x.shape[1] == 602:
            x = self.layer_options['japonica_CK'](x)
        elif x.shape[1] == 635:
            x = self.layer_options['japonica_CD'](x)
        elif x.shape[1] == 85:
            x = self.layer_options['japonica_CE'](x)
        elif x.shape[1] == 333:
            x = self.layer_options['indica_CK'](x)
        elif x.shape[1] == 263:
            x = self.layer_options['indica_CD'](x)
        elif x.shape[1] == 56:
            x = self.layer_options['indica_CE'](x)
        elif x.shape[1] == 743:
            x = self.layer_options['zeamays_CK'](x)
        elif x.shape[1] == 777:
           # print(self.layer_options['zeamays_CD'])
            x = self.layer_options['zeamays_CD'](x)
        elif x.shape[1] == 10:
            x = self.layer_options['zeamays_CE'](x)

        elif x.shape[1] == 1676:
            x = self.layer_options['homo_CK'](x)
        elif x.shape[1] == 1675:
           # print(self.layer_options['zeamays_CD'])
            x = self.layer_options['homo_CD'](x)
        elif x.shape[1] == 151:
            x = self.layer_options['homo_CE'](x)

        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic(space_encoding_tensor)

        rows = edge_index[0]  # 从 edge_index 第一行获取起始节点索引
        cols = edge_index[1]  # 从 edge_index 第二行获取结束节点索引

        # 使用这些行和列索引从 space_encoding_tensor 中提取对应的值
        selected_values = att_bias[rows, cols]

        # 如果您需要保持 selected_values 的形状为 [num_edges, 1] 或类似格式，可以调整形状
        att_bias = selected_values#.unsqueeze(1)
        att_bias = att_bias.unsqueeze(1)
        #ic(att_bias.shape)
        #att_bias = None#self.GraphAttnBias(x,space_encoding_tensor,att_bias)
        #ic(x.shape,adj.shape)
        x = self.GraphNodeFeature(x,adj)
        #ic(x.shape, edge_index.shape, edge_attr.shape)
        x0 = x

        att_bias1 = self.attn_bias1(att_bias)
        #ic(x.shape,edge_index.shape)
        edge_attr1 = self.transform1(x[edge_index[0]] + x[edge_index[1]])
        x1 = F.relu(self.norm1(self.conv1(x, edge_index, edge_attr1,att_bias=att_bias1)))
        x1=x1+self.gate1(x0)*F.sigmoid(x1)
        x = F.dropout(x1, p=0.5)

        # att_bias2 = self.attn_bias2(att_bias)
        # edge_attr2 = self.transform2(x[edge_index[0]] + x[edge_index[1]])
        # x2 = F.relu(self.norm2(self.conv2(x, edge_index, edge_attr2,att_bias=att_bias2)))
        # x2=x2+self.gate2(x)*F.sigmoid(x2)
        # x = F.dropout(x2, p=0.5)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)

        att_bias3 = self.attn_bias3(att_bias)
        edge_attr3 = self.transform3(x[edge_index[0]] + x[edge_index[1]])
        x3 = F.relu(self.norm3(self.conv3(x, edge_index, edge_attr3,att_bias=att_bias3)))
        x3=x3+self.gate3(x)*F.sigmoid(x3)
        # 从节点特征表示中抽取特定训练样本的特征
        #ic(x3.shape)
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]
        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)

        #ic(pred.shape)
        if self.visualize:
            return input,adj,x0,edge_index,pred,tf_embed,target_embed,x1,self.conv1.alpha,x3,self.conv3.alpha
        else:
            return pred
        
from models.GraphormerconvLayer import GraphormerConv
class Graphormer2Link(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, reduction='mean'):
        super(Graphormer2Link, self).__init__()
        dropout=0.1
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.GraphNodeFeature = GraphNodeFeature()
        self.GraphAttnBias = GraphAttnBias(num_heads=3,num_spatial=7+1)
        # 图 Transformer 卷积层
        self.conv1 = GraphormerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv2 = GraphormerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        #self.conv3 = TransformerConv(hidden2_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv3 = GraphormerConv(hidden2_dim*3*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,space_encoding_tensor,att_bias):

        ic(space_encoding_tensor.max())
        att_bias = self.GraphAttnBias(x,space_encoding_tensor,att_bias)
        edge_attr=None
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        x = self.GraphNodeFeature(x,adj)
        x = F.relu(self.conv1(x, edge_index, edge_attr,att_bias=att_bias))
        x = F.dropout(x, p=0.5, training=self.training)
        x = F.relu(self.conv2(x, edge_index, edge_attr,att_bias=att_bias))
        x = F.dropout(x, p=0.5, training=self.training)
        #x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)
        x = self.conv3(x, edge_index, edge_attr,att_bias=att_bias)
        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x[train_sample[:, 0], :]
        target_embed = x[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        return pred
#==================================多种非生物胁迫下的应对机制================================
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import TransformerConv, LayerNorm, global_mean_pool

class StressResponseGRN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim, num_conditions, dropout_rate=0.1):
        super(StressResponseGRN, self).__init__()
        self.num_conditions = num_conditions
        
        # 共享图卷积层
        self.shared_conv1 = TransformerConv(input_dim, hidden_dim, heads=3, dropout=dropout_rate)
        self.shared_conv2 = TransformerConv(hidden_dim * 3, hidden_dim, heads=3, dropout=dropout_rate)
        self.norm1 = LayerNorm(hidden_dim * 3)
        self.norm2 = LayerNorm(hidden_dim * 3)
        
        # 特定于逆境的层
        self.condition_layers = nn.ModuleList([
            nn.Sequential(
                TransformerConv(hidden_dim * 3, hidden_dim, heads=3, dropout=dropout_rate),
                LayerNorm(hidden_dim * 3)
            ) for _ in range(num_conditions)
        ])
        
        # 输出层
        self.output_conv = TransformerConv(hidden_dim * 3, output_dim, heads=1, dropout=dropout_rate, concat=False)
        self.output_norm = LayerNorm(output_dim)

        # 解码器层
        self.decoder = nn.Linear(output_dim, 2)

    def forward(self, x, edge_index, batch_index, condition_index):
        # 共享层
        x = F.relu(self.norm1(self.shared_conv1(x, edge_index)))
        x = F.dropout(x, p=0.1, training=self.training)
        x = F.relu(self.norm2(self.shared_conv2(x, edge_index)))

        # 条件特定层
        condition_outs = []
        for i, condition_layer in enumerate(self.condition_layers):
            condition_mask = (condition_index == i)
            if condition_mask.any():
                condition_data = x[condition_mask]
                condition_edge_index = edge_index[:, condition_mask]
                condition_out = F.relu(condition_layer(condition_data, condition_edge_index))
                condition_outs.append(condition_out)
            else:
                condition_outs.append(None)
        
        # 聚合不同条件的输出
        pooled_outputs = torch.cat([global_mean_pool(cond_out, batch_index[condition_index == i]) 
                                    if cond_out is not None else torch.zeros((1, output_dim)).to(x.device)
                                    for i, cond_out in enumerate(condition_outs)], dim=0)
        
        # 输出层
        output = self.output_norm(self.output_conv(pooled_outputs, edge_index))
        
        # 解码器预测
        predictions = self.decoder(output)
        
        return predictions

# 模型实例化和使用
model = StressResponseGRN(input_dim=100, hidden_dim=64, output_dim=32, num_conditions=4)


from models.graph_transformer_pytorch import GraphTransformer
class Graphormer1Link_alphafold(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, hidden_dim,visualize,reduction,**kwgs):
        super(Graphormer1Link_alphafold, self).__init__()
        dropout=0.01
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize
        adj = kwgs.get('adj',None)
        adj_dense = adj.to_dense()
        out_degree = torch.sum(adj_dense, dim=1).long()
        max_out_degree = torch.max(out_degree)
        in_degree = torch.sum(adj_dense, dim=0).long()
        max_in_degree = torch.max(in_degree)    
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=hidden_dim,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        # 图 Transformer 卷积层

        # self.conv1 = SAGEConv(input_dim, hidden1_dim)
        # self.conv2 = SAGEConv(hidden1_dim, hidden2_dim)
        self.conv1 = GraphTransformer(dim=input_dim, depth=3,edge_dim=edge_dim,with_feedforwards=True,gated_residual=True,rel_pos_emb=True,accept_adjacency_matrix=False)#hidden1_dim, heads=3, dropout=dropout, edge_dim)
        # self.conv2 = GraphTransformer(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = GraphTransformer(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        # self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward(self, x, adj,train_sample,gene_num):
        edge_attr=None
        #图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------
        if self.flag:
            if not adj.is_sparse:
                adj = adj.to_sparse()
            adj_matrix_coalesced = adj.coalesce()
            #adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        #ic('input',x.shape)
        
        x = self.GraphNodeFeature(x,adj)
        #ic('input',x.shape)
        #adj=adj.to_dense()
        #ic('input+edge',x.shape,edge_index.shape)
        mask = torch.ones(1, x.shape[0]).bool()
        x1,_ = F.relu(self.conv1(x[None,], edge_index[None,], mask=mask))
        x1 = torch.squeeze(x1,dim=0)
        #ic('conv1',x1.shape)
        x = F.dropout(x1, p=0.5, training=self.training)
        # x2 = F.relu(self.conv2(x, edge_index, edge_attr))
        # ic('conv2',x2.shape)
        # x = F.dropout(x2, p=0.5, training=self.training)
        # #x = F.relu(self.conv3(x, edge_index, edge_attr))
        # #x = F.dropout(x, p=0.5, training=self.training)
        # x3 = self.conv3(x, edge_index, edge_attr)
        # ic('conv3',x3.shape)
        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x[train_sample[:, 0], :]
        target_embed = x[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize:
            return pred,tf_embed,target_embed,x1,self.conv1.alpha,x2,self.conv2.alpha,x3,self.conv3.alpha
        else:
            return pred

class Graphormer1Link_time(torch.nn.Module):
    def __init__(self, input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type, hidden_dim,visualize,reduction,**kwgs):
        super(Graphormer1Link_time, self).__init__()
        dropout=0.01
        edge_dim=None
        self.flag = flag
        self.device = device
        self.training = training  # 注意：在 PyTorch 中，通常不需要将 training 作为构造函数参数传入，可以使用 model.train() 和 model.eval() 来进行控制
        self.type = type
        self.visualize = visualize
        adj = kwgs.get('adj',None)
        adj_dense = adj.to_dense()
        out_degree = torch.sum(adj_dense, dim=1).long()
        max_out_degree = torch.max(out_degree)
        in_degree = torch.sum(adj_dense, dim=0).long()
        max_in_degree = torch.max(in_degree)    
        self.GraphNodeFeature = GraphNodeFeature(hidden_dim=hidden_dim,num_in_degree=max_in_degree+1,num_out_degree=max_out_degree+1)
        
        # 图 Transformer 卷积层

        # self.conv1 = SAGEConv(input_dim, hidden1_dim)
        # self.conv2 = SAGEConv(hidden1_dim, hidden2_dim)
        self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)

        # self.conv1 = TransformerConv(input_dim, hidden1_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv2 = TransformerConv(hidden1_dim * 3, hidden2_dim, heads=3, dropout=dropout, edge_dim=edge_dim)
        # self.conv3 = TransformerConv(hidden2_dim*3, output_dim, heads=1, dropout=dropout, edge_dim=edge_dim,concat=False)
        self.lstm = nn.LSTM(128,128,batch_first=True)
        self.decode_linear = Linear(output_dim, 2)
        if self.type == 'MLP':
            self.linear = nn.Sequential(
                nn.Linear(2 * output_dim, hidden3_dim),
                nn.ReLU(),
                nn.Linear(hidden3_dim, 1)
            )
            
    def decode(self,tf_embed,target_embed):
        if self.type =='dot':
            prob = torch.mul(tf_embed, target_embed)
            #ic(prob.shape)
            if self.flag:
                prob = torch.sum(prob,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
            else:
                prob = torch.sum(prob,dim=1).view(-1,1)
            #ic(prob.shape)
            return prob
        elif self.type =='cosine':
            if self.flag:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)#self.decode_linear(prob)
                #ic(prob.shape)
                prob_class_0 = 1 - prob
                prob_class_1 = prob
                # 将两个类的概率合并成形状为[256, 2]的张量
                prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                #prob = self.decode_linear(prob)
            else:
                prob = torch.cosine_similarity(tf_embed,target_embed,dim=1).view(-1,1)
            return prob
        elif self.type == 'MLP':
            if self.flag:
                prob = torch.cat([tf_embed, target_embed],dim=1)
                #ic(prob.shape)
                #prob = torch.cat((prob_class_0, prob_class_1), dim=1)
                prob = self.decode_linear(prob)
            else:
                h = torch.cat([tf_embed, target_embed],dim=1)
                prob = self.linear(h)
            #ic(prob.shape)
            return prob
        else:
            raise TypeError(r'{} is not available'.format(self.type))
    
    def forward_1time(self, x, adj):
        edge_attr=None
        # x N , T ,D
        if self.flag:
            if not adj.is_sparse:
                adj = adj.to_sparse()
            adj_matrix_coalesced = adj.coalesce()
            #adj_matrix_coalesced = adj.coalesce()
            # 获取稀疏张量的坐标
            sparse_coords = adj_matrix_coalesced.indices()
            # 有向图的边不需要上三角筛选
            edge_index = sparse_coords
        #----------只适合无向图，邻接矩阵是对称的,取上三角---------------
        else:
            adj_matrix=adj
            adj_matrix_coalesced = adj_matrix.coalesce()
            #ic(adj_matrix_coalesced)
            # 获取聚合稀疏张量的坐标和值
            sparse_coords = adj_matrix_coalesced.indices()
            values = adj_matrix_coalesced.values()
            # 只选择那些上三角部分的坐标（假设邻接矩阵是对称的）
            row_indices, col_indices = sparse_coords[0], sparse_coords[1]
            mask = row_indices < col_indices
            edge_index = torch.stack([row_indices[mask], col_indices[mask]], dim=0)
        # 图 Transformer 需要传入 edge_index 以及 edge_attr（如果边上有特征的话）
        #---------适合有向图，取所有边----------------------------

        ic('input',x.shape)
        x = self.GraphNodeFeature(x,adj)
        #ic('input',x.shape)

                #ic('input+edge',x.shape,edge_index.shape)
        x1 = F.relu(self.conv1(x, edge_index, edge_attr))
        #ic('conv1',x1.shape)
        x = F.dropout(x1, p=0.5, training=self.training)
        x2 = F.relu(self.conv2(x, edge_index, edge_attr))
        ic('conv2',x2.shape)
        x = F.dropout(x2, p=0.5, training=self.training)
        x = F.relu(self.conv3(x, edge_index, edge_attr))
        #x = F.dropout(x, p=0.5, training=self.training)
        x3 = self.conv3(x, edge_index, edge_attr)
        ic('conv3',x3.shape)

        return x3
    
    def forward(self,x, adj,train_sample,gene_num):
       
        N,T,D = x.shape
        xt = [self.forward_1time(x[:,i],adj) for i in range(T)]
        xt = torch.stack(xt,dim=1)
        xt,_ = self.lstm(xt)
        
        x3 = torch.mean(xt,dim=1)
        #x3 = xt[:,-1]

        # 从节点特征表示中抽取特定训练样本的特征
        tf_embed = x3[train_sample[:, 0], :]
        target_embed = x3[train_sample[:, 1], :]

        # 调用解码器根据节点嵌入计算预测
        pred = self.decode(tf_embed, target_embed)
        if self.visualize:
            return pred,tf_embed,target_embed,x1,self.conv1.alpha,x2,self.conv2.alpha,x3,self.conv3.alpha
        else:
            return pred



#---------------------------仅用高光谱数据  nobnCNN------------------------------
#from torchinfo import summary
if __name__ == "__main__":
    input_dim = 30000  # Input feature dimension
    hidden1_dim = 64  # Hidden dimension of the first set of Attention layers
    hidden2_dim = 32  # Hidden dimension of the second set of Attention layers
    hidden3_dim = 10  # Hidden dimension for the tf/target linear networks
    output_dim = 5  # Output feature dimension
    num_head1 = 2  # Number of attention heads in first set
    num_head2 = 3  # Number of attention heads in second set
    alpha = 0.2  # Alpha for the leaky ReLU
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')  # Determine which device to use
    reduction = 'mean'  # Reduction type after attention heads
    type_ = 'dot'  # Decode type for the output layer
    flag = 'True'
    hidden_dim = 30000
    visualize = False
    training = False
    # 创建模型实例并传入必要的参数
    genelink_model =Graphormer1Link(input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim, training,
                 num_head2, flag, device, type_, hidden_dim,visualize).to(device)
    # genelink_model = GENELink(input_dim, hidden1_dim, hidden2_dim, hidden3_dim, output_dim,
    #                         num_head1, num_head2, alpha, device, type_, reduction).to(device)
    # 创建一些假的数据来测试模型
    x = torch.randn((50, input_dim)).to(device)  # 假设有50个样本输入，并且每个样本有 input_dim 个特征
    adj = torch.ones((50, 50)).to(device)  # 假设所有的节点都是互相连接的
    train_sample = torch.randint(high=50, size=(20, 2)).to(device)  # 随机选取20个节点对进行小批量训练
    
    # 前向过程
    pred = genelink_model(x, adj, train_sample)
    
    # 打印模型输出
    # print("Predicted output:", pred)
    # from ptflops import get_model_complexity_info
    # import re
    # macs,params = get_model_complexity_info(net,tuple(x.shape[1:]),as_strings=True,
    #                                         print_per_layer_stat=False,verbose=True)

    # flops = eval(re.findall(r'([\d.]+)',macs)[0])*2
    # # Extract the unit
    # flops_unit = re.findall(r'([A-Za-z]+)',macs)[0][0]

    # print('Computational complexity: {:<8}'.format(macs))
    # print('Computational complexity: {} {}Flops'.format(flops, flops_unit))
    # print('Number of parameters: {:<8}'.format(params))
