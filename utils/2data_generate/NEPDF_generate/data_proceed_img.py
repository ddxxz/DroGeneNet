# Usage: python get_xy_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt gene_pair_list  data_separation_index_list  bulk_expression_data  sc_exprsssion_data flag (0, no label. 1, label)
# command line in developer's linux machine :
# python get_xy_label_data_cnn_combine_from_database.py bulk_gene_list.txt sc_gene_list.txt mmukegg_new_new_unique_rand_labelx_sy.txt mmukegg_new_new_unique_rand_labelx_num_sy.txt /home/yey3/sc_process_1/new_bulk_mouse/prs_calculation/mouse_bulk.h5 /home/yey3/sc_process_1/rank_total_gene_rpkm.h5 1
#################INPUT################################################################################################################################
# 1, bulk_gene_list.txt is the list that convert bulk expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t bulk gene ID'
# 2, sc_gene_list.txt is the list that convert sc expression data gene set into gene symbol IDs. format: 'gene symbol IDs\t sc gene ID'
# 3, gene_pair_list is the list that contains gene pairs and their labels. format : 'GeneA    GeneB     0'
# 4, data_separation index list is a number list that divide gene_pair_list into small parts
# here we use data separation index list to divide gene pairs into small data parts, and make sure that the gene pairs in each index inteval is completely isolated from others. And we can evaluate CNNC's performance on only a small data part.
# if users do not need to separate data, they can just generate a index list to divide the data into N equal parts.
# 5, bulk expression data  it should be a hdf5 format. users can use their own data or data we provided.
# 6, sc expression data  it should be a hdf5 format. users can use their own data or data we provided.
# 7, flag (0, no label. 1, label)
#################OUTPUT
# it generate a data_label folder, and a series of data files containing Nxdata_tf (NEPDF file), ydata_tf (label file) and zdata_tf (gene symbol pair file) for each data part divided.
#python get_xy_label_data_cnn_combine_from_database.py /mnt/h/study/deep_learning/gene/project/CNNC/OsData_FPKM_right/bulk_gene_list.txt None /mnt/h/study/deep_learning/gene/project/CNNC/OsData_FPKM_right/unique_labelx.txt /mnt/h/study/deep_learning/gene/project/CNNC/data/mmukegg_new_new_unique_rand_labelx_num_sy.txt /mnt/h/study/deep_learning/gene/project/CNNC/OsData_FPKM_right/ExpressionData.csv None 1
#python get_xy_label_data_cnn_combine_from_database.py /mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/bulk_gene_list.txt /mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/unique_labelx.txt /mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/labelx_num.txt /mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/ExpressionData.csv 1

import pandas as pd
from numpy import *
import numpy as np
import json, re,os, sys
from icecream import install,ic
install()

save_dir = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/img'#os.path.join(os.getcwd(),'NEPDF_data')
gene_bulk_list = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/bulk_gene_list.txt'
rpkm_bulk_path = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/ExpressionData.csv'
labelx_path = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/unique_labelx.txt'
labelx_num_path = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/labelx_num.txt'
flag = '1'

if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

def get_gene_list_bulk(file_name):
    import re
    h={}
    s = open(file_name,'r')   #gene symbol ID list of bulk RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2)   # h [gene symbol] = gene ID
    s.close()
    return h

def get_gene_list(file_name):
    import re
    h={}
    s = open(file_name,'r') #gene symbol ID list of sc RNA-seq
    for line in s:
        search_result = re.search(r'^([^\s]+)\s+([^\s]+)',line)
        h[search_result.group(1).lower()]=search_result.group(2) # h [gene symbol] = gene ID
    s.close()
    return h

def get_sepration_index(file_name):
    import numpy as np
    index_list = []
    #ic(file_name)
    s = open(file_name, 'r')
    for line in s:
        index_list.append(int(line))
    return (np.array(index_list))

# Script starts from here

h_gene_list_bulk =get_gene_list_bulk(gene_bulk_list) #'bulk_gene_list.txt')#
#print('h_gene_list_bulk',len(h_gene_list_bulk.keys()),list(h_gene_list_bulk.keys())[:3],list(h_gene_list_bulk.values())[:3])
#print ('read bulk gene list')

rpkm_bulk = pd.read_csv(rpkm_bulk_path)#水稻数据

########## generate NEPDF matrix
gene_pair_label = []
s=open(labelx_path)#'mmukegg_new_new_unique_rand_labelx.txt')#)   ### read the gene pair and label file 读取标签
for line in s:
    gene_pair_label.append(line)
#print('gene_pair_label',len(gene_pair_label),gene_pair_label[:3])
gene_pair_index = get_sepration_index(labelx_num_path)#'mmukegg_new_new_unique_rand_labelx_num.npy')#sys.argv[6]) # read file speration index
#print('gene_pair_index',gene_pair_index.shape,gene_pair_index[:3])
s.close()
gene_pair_label_array = array(gene_pair_label)
for i in range(len(gene_pair_index)-1):   #### many sperations
    #print (i)
    start_index = gene_pair_index[i]
    end_index = gene_pair_index[i+1]
    x = []
    y = []
    z = []
    for gene_pair in gene_pair_label_array[start_index:end_index]: ## each speration
        separation = gene_pair.split()
        if flag == '1':
            x_gene_name,y_gene_name,label = separation[0],separation[1],separation[2]
            y.append(label)
        else:
            x_gene_name, y_gene_name = separation[0], separation[1]
        z.append(x_gene_name+'\t'+y_gene_name)


        x_tf_bulk = log10(rpkm_bulk[h_gene_list_bulk[x_gene_name]] + 10 ** -2)  ##[0:249] 249 means the number of samples, users can just remove '[0:249]'
        #ic(y_gene_name)
        #ic(h_gene_list_bulk)
        x_gene_bulk = log10(rpkm_bulk[h_gene_list_bulk[y_gene_name]] + 10 ** -2)#[0:249]
        ic(x_gene_bulk,x_tf_bulk)
        H_T_bulk = histogram2d(x_tf_bulk, x_gene_bulk, bins=32)
        H_bulk= H_T_bulk[0].T
        HT_bulk = (log10(H_bulk / 43261 + 10 ** -4) + 4)/4
        x.append(HT_bulk)
    if (len(x)>0):
        xx = array(x)[:, :, :, newaxis]#增加维度
    else:
        xx = array(x)
    save(save_dir+'/Nxdata_tf' + str(i) + '.npy', xx)
    if flag == '1':
        save(save_dir+'/ydata_tf' + str(i) + '.npy', array(y))
    save(save_dir+'/zdata_tf' + str(i) + '.npy', array(z))



