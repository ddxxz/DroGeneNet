
import matplotlib.pyplot as plt
import os
config = {
    "font.family":'Times New Roman',
    "font.size": 25,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)
def count_csv_files(folder_path, groups,dataset):
    group_counts = {group: 0 for group in groups}
    for group in groups:
        group_folder = os.path.join(folder_path, group,dataset)
        print(group_folder)
        if os.path.exists(group_folder):
            csv_files = [f for f in os.listdir(group_folder) if f.endswith('.csv')]
            group_counts[group] = len(csv_files)
    print(group_counts)
    return group_counts

# 示例文件夹路径和组别

data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type== 'KEGG':
        groups = [f'{data_type}_indica',f'{data_type}_japonica' ,f'{data_type}_zeamays' ]
    else:
        groups = [f'{data_type}_indica_BGI',  f'{data_type}_japonica',f'{data_type}_zeamays']

    folder_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/louvein_GRNcluster/'
    # 获取两个数据集的CSV文件数量
    counts1 = count_csv_files(folder_path , groups,'CK-CD')
    counts2 = count_csv_files(folder_path, groups,'CK-CE')

    # 绘制图表
    fig, axs = plt.subplots(1, 2, figsize=(10, 5))
    custom_labels = ['indica', 'japonica', 'zeamays']
    colors1 = ['#1f77b4', '#1f77b4', '#1f77b4']  # 这里可以放置你喜欢的颜色
    colors2 = ['#8c564b', '#8c564b', '#8c564b']
    # 第一个子图
    axs[0].bar(groups, [counts1[group] for group in groups],color=colors1)
    axs[0].set_title(f'{data_type} CK CD', fontsize=35, fontweight='bold')
    #axs[0].set_xlabel('Groups')
    axs[0].set_ylabel('Number of module', fontsize=30, fontweight='bold')
    axs[0].set_xticks(range(len(groups)))
    axs[0].set_xticklabels(custom_labels)
    plt.setp(axs[0].get_xticklabels(), rotation=30, ha='right')
    # 第二个子图
    axs[1].bar(groups, [counts2[group] for group in groups],color=colors2)
    axs[1].set_title(f'{data_type} CK CE', fontsize=35, fontweight='bold')
    #axs[1].set_xlabel('Groups')
    axs[1].set_ylabel('Number of module', fontsize=30, fontweight='bold')
    axs[1].set_xticks(range(len(groups)))
    axs[1].set_xticklabels(custom_labels)
    plt.setp(axs[1].get_xticklabels(), rotation=30, ha='right')
    plt.tight_layout()

    # 显示图表
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/louvein_GRNcluster/{data_type}_modulenum.png',bbox_inches='tight')







