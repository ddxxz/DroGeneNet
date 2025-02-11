import pandas as pd
from icecream import install,ic
install()
import seaborn as sns
import matplotlib.pyplot as plt
config = {
    "font.family":'Times New Roman',
    "font.size": 22,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)

def remove_extremes(group):
    # 排除最高和最低值
    return group[(group['pred_weights2'] != group['pred_weights2'].max()) & 
                 (group['pred_weights2'] != group['pred_weights2'].min())]
titles=['indica','zeamays','japonica',]
data_types=['PPI']#'TF',,'KEGG'
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']
    datapath='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight'
    for i,speci in enumerate(species):
        CK_data=pd.read_csv(f'{datapath}/{data_type}_{speci}/CK_all.txt',sep='\t')
        CD_data=pd.read_csv(f'{datapath}/{data_type}_{speci}/CD_all.txt',sep='\t')
        ic(CK_data.keys())
        CK_TF = CK_data.groupby('tf_family')['pred_weights2']
        CD_TF = CD_data.groupby('tf_family')['pred_weights2']
        ic(CK_TF)
        CK_data['Condition'] = 'CK'
        CD_data['Condition'] = 'CD'

        CK_data = CK_data.groupby('tf_family').apply(remove_extremes).reset_index(drop=True)
        CD_data = CD_data.groupby('tf_family').apply(remove_extremes).reset_index(drop=True)

        CK_data = CK_data[['tf_family', 'pred_weights2', 'Condition']]
        CD_data = CD_data[['tf_family', 'pred_weights2', 'Condition']]
        
        # 合并CK和CD数据
        combined_data = pd.concat([CK_data, CD_data])
        colors = {"CK": "#00798c", "CD": "#d1495b"}
        # 为每个转录因子家族绘制密度图
        for family in combined_data['tf_family'].unique():
            plt.figure(figsize=(5, 5))
            subset = combined_data[combined_data['tf_family'] == family]
            #ax = sns.violinplot(data=subset, x='Condition', y='pred_weights2', palette=colors)
            sns.boxplot(data=subset, x='Condition', y='pred_weights2', palette=colors)
            # ax=sns.kdeplot(data=subset, x='pred_weights2', hue='Condition', fill=True, common_norm=False, alpha=0.5, palette=colors)


            # 添加图标题和坐标轴标签
            plt.title(f'{titles[i]} {family}',fontsize=30, weight='bold', pad=15)
            plt.xlabel('Group',fontsize=25, weight='bold')
            if data_type== 'KEGG' or data_type== 'PPI':
                plt.ylabel('protein interaction activity',fontsize=25, weight='bold')
            else:
                plt.ylabel('TF activity',fontsize=25, weight='bold')
            plt.gca().spines['top'].set_linewidth(1.5)
            plt.gca().spines['bottom'].set_linewidth(1.5)
            plt.gca().spines['right'].set_linewidth(1.5)
            plt.gca().spines['left'].set_linewidth(1.5)
            # 调整图例位置和大小
            #handles, labels = ax.get_legend_handles_labels()

            # 使用获取的句柄和标签手动设置图例
            #ax.legend(handles=handles, labels=labels, title='Condition', loc='upper left', fontsize=16, title_fontsize=18)

            #plt.legend(title=f'{family} TF Family', bbox_to_anchor=(1.05, 1), loc=2)
            plt.tight_layout()
            if '/' in family:
                family=family.split('/')[0]
            plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_activity/{data_type}_{speci}_{family}_CK_CD.png',bbox_inches='tight')
            plt.close()


    for i,speci in enumerate(species):
        CK_data=pd.read_csv(f'{datapath}/{data_type}_{speci}/CK_all.txt',sep='\t')
        CD_data=pd.read_csv(f'{datapath}/{data_type}_{speci}/CD_all.txt',sep='\t')
        CE_data=pd.read_csv(f'{datapath}/{data_type}_{speci}/CE_all.txt',sep='\t')
        ic(CK_data.keys())
        CK_TF = CK_data.groupby('tf_family')['pred_weights2']
        CD_TF = CD_data.groupby('tf_family')['pred_weights2']
        CE_TF = CE_data.groupby('tf_family')['pred_weights2']
        ic(CK_TF)
        CK_data['Condition'] = 'CK'
        CD_data['Condition'] = 'CD'
        CE_data['Condition'] = 'CE'

        CK_data = CK_data.groupby('tf_family').apply(remove_extremes).reset_index(drop=True)
        CD_data = CD_data.groupby('tf_family').apply(remove_extremes).reset_index(drop=True)
        CE_data = CE_data.groupby('tf_family').apply(remove_extremes).reset_index(drop=True)

        CK_data = CK_data[['tf_family', 'pred_weights2', 'Condition']]
        CD_data = CD_data[['tf_family', 'pred_weights2', 'Condition']]
        CE_data = CE_data[['tf_family', 'pred_weights2', 'Condition']]
        
        # 合并CK和CD数据
        combined_data = pd.concat([CK_data, CD_data,CE_data])
        colors = {"CK": "#00798c", "CD": "#d1495b","CE": "#edae49"}
        # 为每个转录因子家族绘制密度图
        for family in combined_data['tf_family'].unique():
            plt.figure(figsize=(5, 5))
            subset = combined_data[combined_data['tf_family'] == family]
            #ax = sns.violinplot(data=subset, x='Condition', y='pred_weights2', palette=colors)
            sns.boxplot(data=subset, x='Condition', y='pred_weights2', palette=colors)
            # ax=sns.kdeplot(data=subset, x='pred_weights2', hue='Condition', fill=True, common_norm=False, alpha=0.5, palette=colors)


            # 添加图标题和坐标轴标签
            plt.title(f'{titles[i]} {family}',fontsize=30, weight='bold', pad=15)
            plt.xlabel('Group',fontsize=25, weight='bold')
            if data_type== 'KEGG' or data_type== 'PPI':
                plt.ylabel('protein interaction activity',fontsize=25, weight='bold')
            else:
                plt.ylabel('TF activity',fontsize=25, weight='bold')
            plt.gca().spines['top'].set_linewidth(1.5)
            plt.gca().spines['bottom'].set_linewidth(1.5)
            plt.gca().spines['right'].set_linewidth(1.5)
            plt.gca().spines['left'].set_linewidth(1.5)
            # 调整图例位置和大小
            #handles, labels = ax.get_legend_handles_labels()

            # 使用获取的句柄和标签手动设置图例
            #ax.legend(handles=handles, labels=labels, title='Condition', loc='upper left', fontsize=16, title_fontsize=18)

            #plt.legend(title=f'{family} TF Family', bbox_to_anchor=(1.05, 1), loc=2)
            plt.tight_layout()
            if '/' in family:
                family=family.split('/')[0]
            plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_activity/{data_type}_{speci}_{family}_CK_CD_CE.png',bbox_inches='tight')
            plt.close()














































