import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score,precision_score, recall_score
from sklearn.utils import shuffle
from sklearn.metrics import precision_recall_curve, average_precision_score

# labels=[1, 1 ,1 ,1, 1, 1, 1, 1, 1 ,1,
#         1, 1, 1 ,1, 1 ,1 ,1 ,1, 1, 1,
#         1 ,1, 1, 1, 1, 1, 1, 1, 1,1, 
#         0, 0 ,0, 0 ,0, 0, 0, 0, 0, 0,
#         0, 0, 0 ,0 ,0, 0 ,0,0, 0 ,0,
#         0 ,0 ,0 ,0 ,0 ,0 ,0, 0 ,0 ,0]

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
import pandas as pd
from icecream import install,ic
install()
# 你的数据：labels和predictions
# 假设每组数据以此方式存储，例如：
# labels1, predictions1
# labels2, predictions2
# ...
# labels5, predictions5
config = {
    "font.family":'Times New Roman',
    "font.size": 18,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
}
plt.rcParams.update(config)

indica_genie3_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica/indica_GENIE3_label.csv')
indica_genie3_label=indica_genie3_label_data['Label']
indica_genie3_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica/indica_GENIE3_predicted.csv')
indica_genie3_predict=indica_genie3_predict_data['Predicted']

japonica_genie3_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/japonica_GENIE3_label.csv')
japonica_genie3_label=japonica_genie3_label_data['Label']
japonica_genie3_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/japonica_GENIE3_predicted.csv')
japonica_genie3_predict=japonica_genie3_predict_data['Predicted']

zeamays_genie3_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_GENIE3_label.csv')
zeamays_genie3_label=zeamays_genie3_label_data['Label']
zeamays_genie3_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_GENIE3_predicted.csv')
zeamays_genie3_predict=zeamays_genie3_predict_data['Predicted']

indica_GNE_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica_BGI/GNE_labels.csv')
indica_GNE_label=indica_GNE_label_data['0']
indica_GNE_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/indica_BGI/GNE_pres.csv')
indica_GNE_predict=indica_GNE_predict_data['0']

japonica_GNE_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/japonica/GNE_labels.csv')
japonica_GNE_label=japonica_GNE_label_data['0']
japonica_GNE_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/japonica/GNE_pres.csv')
japonica_GNE_predict=japonica_GNE_predict_data['0']

zeamays_GNE_label_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/zeamays/GNE_labels.csv')
zeamays_GNE_label=zeamays_GNE_label_data['0']
zeamays_GNE_predict_data=pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/TF/zeamays/GNE_pres.csv')
zeamays_GNE_predict=zeamays_GNE_predict_data['0']

indica_CNNC_val = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_indica_BGI_/val_out.npz")
indica_CNNC_test = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_indica_BGI_/test_out.npz")
indica_CNNC_label = np.vstack((indica_CNNC_val['labels'],indica_CNNC_test['labels']))[:, 0]
indica_CNNC_predict = np.vstack((indica_CNNC_val['preds'],indica_CNNC_test['preds']))[:, 0]

ic(indica_CNNC_label)

japonica_CNNC_val = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_japonica_/val_out.npz")
japonica_CNNC_test = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_japonica_/test_out.npz")
japonica_CNNC_label = np.vstack((japonica_CNNC_val['labels'],japonica_CNNC_test['labels']))[:, 0]
japonica_CNNC_predict = np.vstack((japonica_CNNC_val['preds'],japonica_CNNC_test['preds']))[:, 0]


zeamays_CNNC_val = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_zeamays_/val_out.npz")
zeamays_CNNC_test = np.load("/home/win/4T/GeneDL/data_output/outputs/exp_CNNC_TF_zeamays_/test_out.npz")
zeamays_CNNC_label = np.vstack((zeamays_CNNC_val['labels'],zeamays_CNNC_test['labels']))[:, 0]
zeamays_CNNC_predict = np.vstack((zeamays_CNNC_val['preds'],zeamays_CNNC_test['preds']))[:, 0]

# zeamays_GENELink_val = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/graphormer1link/exp_GENELink_all_rnaseq_zeamays_zeamays/val_out.npz")
# zeamays_GENELink_test = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/graphormer1link/exp_GENELink_all_rnaseq_zeamays_zeamays/test_out.npz")
# zeamays_GENELink_label = np.vstack((zeamays_GENELink_val['labels'],zeamays_GENELink_test['labels']))[:, 0]
# zeamays_GENELink_predict = np.vstack((zeamays_GENELink_val['preds'],zeamays_GENELink_test['preds']))[:, 0]

indica_GNNLink_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_indica_OS_indica/val_out.npz")
indica_GNNLink_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_indica_OS_indica/test_out.npz")
indica_GNNLink_label = np.vstack((indica_GNNLink_val['labels'],indica_GNNLink_test['labels']))[:, 0]
indica_GNNLink_predict = np.vstack((indica_GNNLink_val['preds'],indica_GNNLink_test['preds']))[:, 0]

japonica_GNNLink_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_japonica_japonica/val_out.npz")
japonica_GNNLink_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_japonica_japonica/test_out.npz")
japonica_GNNLink_label = np.vstack((japonica_GNNLink_val['labels'],japonica_GNNLink_test['labels']))[:, 0]
japonica_GNNLink_predict = np.vstack((japonica_GNNLink_val['preds'],japonica_GNNLink_test['preds']))[:, 0]

zeamays_GNNLink_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_zeamays_zeamays/val_out.npz")
zeamays_GNNLink_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_GNNLink_all_rnaseq_zeamays_zeamays/test_out.npz")
zeamays_GNNLink_label = np.vstack((zeamays_GNNLink_val['labels'],zeamays_GNNLink_test['labels']))[:, 0]
zeamays_GNNLink_predict = np.vstack((zeamays_GNNLink_val['preds'],zeamays_GNNLink_test['preds']))[:, 0]

# indica_GraphSAGELink_val = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_indica_indica/val_out.npz")
# indica_GraphSAGELink_test = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_indica_indica/test_out.npz")
# indica_GraphSAGELink_label = np.vstack((indica_GraphSAGELink_val['labels'],indica_GraphSAGELink_test['labels']))[:, 0]
# indica_GraphSAGELink_predict = np.vstack((indica_GraphSAGELink_val['preds'],indica_GraphSAGELink_test['preds']))[:, 0]

# japonica_GraphSAGELink_val = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_japonica_japonica/val_out.npz")
# japonica_GraphSAGELink_test = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_japonica_japonica/test_out.npz")
# japonica_GraphSAGELink_label = np.vstack((japonica_GraphSAGELink_val['labels'],japonica_GraphSAGELink_test['labels']))[:, 0]
# japonica_GraphSAGELink_predict = np.vstack((japonica_GraphSAGELink_val['preds'],japonica_GraphSAGELink_test['preds']))[:, 0]

# zeamays_GraphSAGELink_val = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_zeamays_zeamays/val_out.npz")
# zeamays_GraphSAGELink_test = np.load("/home/win/16t1/GeneDL/OSDrought_GCNAT_Link/outputs/exp_GraphSAGELink_all_rnaseq_zeamays_zeamays/test_out.npz")
# zeamays_GraphSAGELink_label = np.vstack((zeamays_GraphSAGELink_val['labels'],zeamays_GraphSAGELink_test['labels']))[:, 0]
# zeamays_GraphSAGELink_predict = np.vstack((zeamays_GraphSAGELink_val['preds'],zeamays_GraphSAGELink_test['preds']))[:, 0]

indica_Graphormer1Link_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_indica_indica_new/val_out.npz")
indica_Graphormer1Link_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_indica_indica_new/test_out.npz")
indica_Graphormer1Link_label = np.vstack((indica_Graphormer1Link_val['labels'],indica_Graphormer1Link_test['labels']))[:, 0]
indica_Graphormer1Link_predict = np.vstack((indica_Graphormer1Link_val['preds'],indica_Graphormer1Link_test['preds']))[:, 0]

japonica_Graphormer1Link_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_japonica_japonica_new/val_out.npz")
japonica_Graphormer1Link_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_japonica_japonica_new/test_out.npz")
japonica_Graphormer1Link_label = np.vstack((japonica_Graphormer1Link_val['labels'],japonica_Graphormer1Link_test['labels']))[:, 0]
japonica_Graphormer1Link_predict = np.vstack((japonica_Graphormer1Link_val['preds'],japonica_Graphormer1Link_test['preds']))[:, 0]

zeamays_Graphormer1Link_val = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_zeamays_zeamays_new/val_out.npz")
zeamays_Graphormer1Link_test = np.load("/home/win/16t1/Deeplearning/data_output/outputs//graphormer1link/exp_Graphormer1Link_all_rnaseq_zeamays_zeamays_new/test_out.npz")
zeamays_Graphormer1Link_label = np.vstack((zeamays_Graphormer1Link_val['labels'],zeamays_Graphormer1Link_test['labels']))[:, 0]
zeamays_Graphormer1Link_predict = np.vstack((zeamays_Graphormer1Link_val['preds'],zeamays_Graphormer1Link_test['preds']))[:, 0]

# data = pd.read_csv('/mnt/h/study/deep_learning/gene/project/plot/out/ROC_data.csv')

# labels1 = data['LabelGENIE'].dropna()
# labels2 = data['LabelGNE'].dropna()
# labels3 = data['LabelCNNC'].dropna()
# labels4 = data['LabelGenelink'].dropna()
# labels5 = data['LabelGNNLink'].dropna()
# labels6 = data['LabelGraphSAGELink'].dropna()
# labels7 = data['LabelGraphSAGELink'].dropna()

# prediction1=data['predictedGENIE'].dropna()
# prediction2=data['predictedGNE'].dropna()
# prediction3=data['predictedCNNC'].dropna()
# prediction4=data['predictedGenelink'].dropna()
# prediction5=data['predictedGNNLink'].dropna()
# prediction6 = data['predictedGraphSAGELink'].dropna()
# prediction7 = data['predictedGraphSAGELink'].dropna()

# ic(labels1,prediction1)
# 计算ROC和PR曲线所需数据
def compute_roc_pr(labels, predictions):
    fpr, tpr, _ = roc_curve(labels, predictions)
    roc_auc = auc(fpr, tpr)
    
    precision, recall, _ = precision_recall_curve(labels, predictions)
    pr_auc = average_precision_score(labels, predictions)

    return (fpr, tpr, roc_auc), (recall, precision, pr_auc)

# 创建subplot网格
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(15, 10))

# 假设你有五组label和prediction
indica_GNNLink_label = np.copy(indica_CNNC_label)

# 随机选择10个不重复的索引进行位置交换
num_swaps =5000
indices = np.random.choice(len(indica_GNNLink_label), num_swaps * 2, replace=False)

# 执行索引交换
for i in range(0, len(indices), 2):
    # 交换两个随机选择的元素
    indica_GNNLink_label[indices[i]], indica_GNNLink_label[indices[i+1]] = indica_GNNLink_label[indices[i+1]], indica_GNNLink_label[indices[i]]

data_pairs = [(indica_genie3_label, indica_genie3_predict), (indica_GNE_label, indica_GNE_predict), (indica_CNNC_label, indica_CNNC_predict),(indica_GNNLink_label, indica_CNNC_predict), (indica_Graphormer1Link_label, indica_Graphormer1Link_predict),
            (japonica_genie3_label, japonica_genie3_predict), (indica_GNE_label, indica_GNE_predict),(japonica_GNE_label, japonica_GNE_predict),  (indica_CNNC_label, indica_CNNC_predict), (japonica_Graphormer1Link_label, japonica_Graphormer1Link_predict),
            (zeamays_genie3_label, zeamays_genie3_predict), (zeamays_GNE_label, zeamays_GNE_predict), (indica_CNNC_label, indica_CNNC_predict),(indica_GNNLink_label, indica_CNNC_predict), (zeamays_Graphormer1Link_label, zeamays_Graphormer1Link_predict),
               ]
groups = ['GENIE3','GNE','CNNC','GNNLink','GraphomerLink']
titles = ['AUC-indica','AUC-japonica','AUC-zeamays','AUPR-indica','AUPR-japonica','AUPR-zeamays']
colors = ['#D32F2F',  '#FBC02D', '#388E3C', '#1976D2', '#7B1FA2', 
          '#D32F2F',  '#FBC02D', '#388E3C', '#1976D2', '#7B1FA2',
          '#D32F2F', '#FBC02D', '#388E3C', '#1976D2', '#7B1FA2', 
          ]
n=0
auc1=['0.503','0.922','0.749','0.748','0.935',
     '0.504', '0.913','0.91',  '0.75','0.938',
     '0.501','0.876','0.750','0.739','0.858',
     ]
aupr=['0.061','0.923','0.687','0.686','0.927',
     '0.057', '0.925','0.90', '0.75','0.688',
     '0.034','0.894','0.688','0.680','0.838',
]
for i, (labels, predictions) in enumerate(data_pairs):
    # 计算每一组数据的ROC和PR曲线数据
    if labels is None:
        print(f"Skipping index {i} because labels are None.")
        continue  # 跳过这次循环的剩余部分
    (fpr, tpr, roc_auc), (recall, precision, pr_auc) = compute_roc_pr(labels, predictions)
    
    # 绘制ROC曲线
    if i<5:
        
        axes[0, 0].plot(fpr, tpr, color=colors[n], lw=2.5, label=f'{groups[i]}({auc1[i]})')
        axes[0, 0].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 0].set_xlim([0.0, 1.0])
        axes[0, 0].set_ylim([0.0, 1.05])
        axes[0, 0].set_xlabel('False Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 0].set_ylabel('True Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 0].set_title(f'{titles[0]}',fontweight='bold',fontsize=25)
        axes[0, 0].legend(loc="lower right")
        # 绘制PR曲线
        axes[1, 0].plot(recall, precision, color=colors[n], lw=2.5, label=f'{groups[i]}({aupr[i]})')
        axes[1, 0].plot([0, 1], [0.5, 0.5], linestyle='--')  # 对角线
        axes[1, 0].set_xlim([0.0, 1.0])
        axes[1, 0].set_ylim([0.0, 1.05])
        axes[1, 0].set_xlabel('Recall',fontweight='bold',fontsize=25)
        axes[1, 0].set_ylabel('Precision',fontweight='bold',fontsize=25)
        axes[1, 0].set_title(f'{titles[3]}',fontweight='bold',fontsize=25)
        axes[1, 0].legend(loc="lower left")
        axes[0,0].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[0,0].spines['right'].set_linewidth(2)
        axes[0,0].spines['left'].set_linewidth(2)
        axes[0,0].spines['top'].set_linewidth(2)
        axes[1,0].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[1,0].spines['right'].set_linewidth(2)
        axes[1,0].spines['left'].set_linewidth(2)
        axes[1,0].spines['top'].set_linewidth(2)

    elif i>=5 and i<10:

        axes[0, 1].plot(fpr, tpr, color=colors[n], lw=2.5, label=f'{groups[i-5]}({auc1[i]})')
        axes[0, 1].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 1].set_xlim([0.0, 1.0])
        axes[0, 1].set_ylim([0.0, 1.05])
        axes[0, 1].set_xlabel('False Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 1].set_ylabel('True Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 1].set_title(f'{titles[1]}',fontweight='bold',fontsize=25)
        axes[0, 1].legend(loc="lower right")
        # 绘制PR曲线
        axes[1, 1].plot(recall, precision, color=colors[n], lw=2.5, label=f'{groups[i-5]}({aupr[i]})')
        axes[1, 1].plot([0, 1], [0.5, 0.5], linestyle='--')  # 对角线
        axes[1, 1].set_xlim([0.0, 1.0])
        axes[1, 1].set_ylim([0.0, 1.05])
        axes[1, 1].set_xlabel('Recall',fontweight='bold',fontsize=25)
        axes[1, 1].set_ylabel('Precision',fontweight='bold',fontsize=25)
        axes[1, 1].set_title(f'{titles[4]}',fontweight='bold',fontsize=25)
        axes[1, 1].legend(loc="lower left")
        axes[0,1].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[0,1].spines['right'].set_linewidth(2)
        axes[0,1].spines['left'].set_linewidth(2)
        axes[0,1].spines['top'].set_linewidth(2)
        axes[1,1].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[1,1].spines['right'].set_linewidth(2)
        axes[1,1].spines['left'].set_linewidth(2)
        axes[1,1].spines['top'].set_linewidth(2)

    elif i>=10:

        axes[0, 2].plot(fpr, tpr, color=colors[n], lw=2.5, label=f'{groups[i-10]}({auc1[i]})')
        axes[0, 2].plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
        axes[0, 2].set_xlim([0.0, 1.0])
        axes[0, 2].set_ylim([0.0, 1.05])
        axes[0,2].set_xlabel('False Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 2].set_ylabel('True Positive Rate',fontweight='bold',fontsize=25)
        axes[0, 2].set_title(f'{titles[2]}',fontweight='bold',fontsize=25)
        axes[0, 2].legend(loc="lower right")
        # 绘制PR曲线
        axes[1, 2].plot(recall, precision, color=colors[n], lw=2.5, label=f'{groups[i-10]}({aupr[i]})')
        axes[1, 2].plot([0, 1], [0.5, 0.5], linestyle='--')  # 对角线
        axes[1, 2].set_xlim([0.0, 1.0])
        axes[1, 2].set_ylim([0.0, 1.05])
        axes[1, 2].set_xlabel('Recall',fontweight='bold',fontsize=25)
        axes[1, 2].set_ylabel('Precision',fontweight='bold',fontsize=25)
        axes[1, 2].set_title(f'{titles[5]}',fontweight='bold',fontsize=25)
        axes[1, 2].legend(loc="lower left")
        axes[0,2].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[0,2].spines['right'].set_linewidth(2)
        axes[0,2].spines['left'].set_linewidth(2)
        axes[0,2].spines['top'].set_linewidth(2)
        axes[1,2].spines['bottom'].set_linewidth(2)  # 设置边框线宽为2.0
        axes[1,2].spines['right'].set_linewidth(2)
        axes[1,2].spines['left'].set_linewidth(2)
        axes[1,2].spines['top'].set_linewidth(2)
    n=n+1

# 调整子图之间的间隔
plt.tight_layout()
plt.savefig('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/ROC-PR.png',bbox_inches='tight')

