import omicverse as ov
ov.utils.ov_plot_set()
import matplotlib.pyplot as plt
import pandas as pd
import pickle
import numpy as np
from icecream import install,ic
install()
import pandas as pd
import numpy as np
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import seaborn as sns
#https://www.cnblogs.com/wangshicheng/p/11121537.html
#https://www.jianshu.com/p/901e7dbbcf7d
def has_nan(array):
    return np.isnan(array).any()

#--------------------------------WGCNA分为表达量聚类分析和表型关联两部分------------------------------------

all_species=['indica_BGI','japonica','zeamays']  #
for species in all_species:

    with open(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/{species}_label_sample_expression.pkl','rb') as fp:
        data=pickle.load(fp)

    '''示例数据格式 标准化之后的数据
    rehydration1,rehydration2,dyhydration1,dyhydration2,control1,control2
    Os01g0100100,19.0,35.0,9.0,7.0,35.0,11.0
    Os01g0100100,19.0,35.0,9.0,7.0,35.0,11.0
    Os01g0100100,19.0,35.0,9.0,7.0,35.0,11.0
    Os01g0100100,19.0,35.0,9.0,7.0,35.0,11.0
    Os01g0100100,19.0,35.0,9.0,7.0,35.0,11.0
    '''

    data.columns = data.iloc[0]  # 将第一行数据设置为列名
    # 删除原来的第一行数据
    data = data[1:]
    # 重置索引
    #data.reset_index(drop=True, inplace=True)
    ic(data)
    #data.set_index(0, inplace=True)#FEATURE_ID
    

    df_ck =pd.DataFrame(data['CK'],index=data.index)
    df_cd =pd.DataFrame(data['CD'],index=data.index)
    #df_ce =pd.DataFrame(data['CE'],index=data.index)
    #df_ck.to_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_WGCNA/{species}_CK.csv',sep='\t')
    #df_cd.to_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_WGCNA/{species}_CD.csv',sep='\t')
    ic(df_ck)

    treat=['CK','CD'] #,'CE'

    for i,matrix in enumerate([df_ck,df_cd]):  #,df_ce
        ic(matrix)
        matrix.fillna(0, inplace=True)
        # 删除重复行
        matrix.drop_duplicates(inplace=True)
        matrix.to_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_WGCNA/{species}_{treat[i]}_WGCNA.txt',sep='\t')
        
    #     matrix += 1e-16
    #     ic(matrix)


        # gene_wgcna=ov.bulk.pyWGCNA(matrix,save_path=f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_WGCNA/{species}')
        # #gene_wgcna.calculate_correlation_direct(method='pearson',save=False)   #直接相关性网络
        # gene_wgcna.result = pd.DataFrame(np.corrcoef(gene_wgcna.data.T.to_numpy(),rowvar=False,dtype='float32'),
        #                                 columns=matrix.T.columns,index=matrix.T.columns
        #                                 )
        # #print('result',has_nan(gene_wgcna.result))
        # #将直接相关矩阵转换为间接相关矩阵来计算软阈值，软阈值可以帮助我们将原来的相关网络转换为无尺度网络 
        # #无尺度网络的典型特征是在网络中的大部分节点只和很少节点连接，而有极少的节点与非常多的节点连接
        # gene_wgcna.calculate_correlation_indirect(save=False)
        # #print('temp',has_nan(gene_wgcna.temp))
        # gene_wgcna.calculate_soft_threshold(save=False)
        # #R2越接近1，网络就越接近无尺度网络  利用 β 值，我们可以根据方程将相关矩阵转换成邻接矩阵

        # plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/expression_WGCNA/{species}_{treat[i]}_WGCNA_threshold.png',bbox_inches='tight')
        # #构造拓扑重叠矩阵   !!!注意此处修改了源代码
        # adj = gene_wgcna.calculate_corr_matrix()
        # ic(adj)
        #print('cor',has_nan(gene_wgcna.cor))
        #在获得基因间的拓扑重叠矩阵后，我们使用动态剪切树的方式来寻找基因间的模块。在这里，我们使用WGCNA作者发表的DynamicTree的算法来实现
        # gene_wgcna.calculate_distance()
        # #print('distances',has_nan(gene_wgcna.distances))
        # gene_wgcna.calculate_geneTree()
        # gene_wgcna.calculate_dynamicMods()
        # module=gene_wgcna.calculate_gene_module()
        # plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/expression_WGCNA/{species}_{treat[i]}_WGCNA_co-expressions_module.png',bbox_inches='tight')
        #print(module.head())
        #==========================================================================================
        #可视化拓扑重叠矩阵与模块之间的关系
        # gene_wgcna.plot_matrix()
        # plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/expression_WGCNA/{species}_{treat[i]}_WGCNA.png',bbox_inches='tight')
        # #print(gene_wgcna.get_sub_module([6,12]).shape)
        # #对一个基因或一个通路的模块感兴趣，我们需要提取基因的子模块进行分析和定位
        # gene_wgcna.get_sub_network([6,12],correlation_threshold=0.95)
        # #可视化基因相关性网络
        # gene_wgcna.plot_sub_network([6,12],pos_type='kamada_kawai',pos_scale=3,pos_dim=2,
        #                         figsize=(8,8),node_size=20,label_fontsize=8,
        #                         label_bbox={"ec": "white", "fc": "white", "alpha": 0.6})
        # plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/expression_WGCNA/{species}_{treat[i]}_WGCNA_Sub_co-expression_module.png',bbox_inches='tight')

#-------------------------------------------表型相关性分析(模块和性状相关性分析)----------------------------------------
# meta=ov.utils.read_csv(filepath_or_buffer='/mnt/g/gene/course_try/Omicverse/genesets/character.csv',
#                           index_col=0)
# print(meta.head())
# cor_matrix=gene_wgcna.analysis_meta_correlation(meta)
# ax=gene_wgcna.plot_meta_correlation(cor_matrix)
# plt.savefig('/mnt/g/gene/course_try/Omicverse/out/WGCNA_correlation.png',bbox_inches='tight')

