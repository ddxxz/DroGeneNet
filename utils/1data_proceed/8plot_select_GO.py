# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# from icecream import install,ic
# install()
# config = {
#     "font.family":'Times New Roman',
#     "font.size": 25,
#     "mathtext.fontset":'stix',
#     "font.serif": ['Times New Roman'],
#     "font.weight": 'bold'
# }
# plt.rcParams.update(config)
# # 读取数据
# bp_data = pd.read_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/homo_BP_GO.csv')
# cc_data = pd.read_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/homo_CC_GO.csv')
# mf_data = pd.read_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/homo_MF_GO.csv')

# combined_df = pd.concat([bp_data, cc_data, mf_data], ignore_index=True)

# # 定义处理和绘图的函数
# def process_and_plot(data, title):
#     # 筛选Adjusted P-value小于0.05的行
#     filtered_data = data[data['Adjusted P-value'] < 0.05]
#     #filtered_data.to_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/homo_select_GO.csv',index=False)
#     filtered_data = pd.read_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/homo_select_GO.csv')
#     # 设置横坐标为Overlap中的分数，这里我们将它转换为整数分子以供后续使用
#     richfactor = filtered_data['Overlap']
#     richfactor_percent = richfactor.apply(lambda x: float(x.split('/')[0]) / float(x.split('/')[1]) * 100 if '/' in x else None)
#     numerators = richfactor.apply(lambda x: int(x.split('/')[0]) if '/' in x and x.split('/')[0].isdigit() else None)

#     # 设置基因数量表示气泡大小
#     sizes = numerators

#     # y轴设置为Term列的条目
#     terms = filtered_data['Term']

#     # 设置颜色映射Adjusted P-value
#     colors = filtered_data['Adjusted P-value']

#     # 绘图
#     fig, ax = plt.subplots(figsize=(5, 7))
#     scatter = ax.scatter(richfactor_percent, terms, s=sizes * 300, c=colors, cmap='RdBu', edgecolors='black', linewidth=1, vmin=0, vmax=0.05,zorder=3) #alpha=0.6, 
#     ax.set_title(title, fontsize=35, fontweight='bold',pad=10)
#     ax.set_xlabel('Rich Factor (%)',fontsize=30,fontweight='bold') 
#     ax.set_ylabel('GO Term',fontsize=30,fontweight='bold')
#     ax.set_xticks([0, 25, 50, 75, 100])  # 设置x轴的固定刻度
#     ax.set_xlim(-10, 110)
#     #ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
#     # 绘制横向网格线
#     ax.grid(True, axis='y', linestyle='--', linewidth=0.7, zorder=1)
#     # 调整颜色条的位置，使其上下对半分

#     cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
#     #cbar.set_label('FDR', rotation=0, labelpad=10)
#     fig.text(0.77, 0.45, 'Q value', rotation=0, ha='left', va='center', fontsize=22)
#     cbar.ax.xaxis.set_label_position('top')  # 设置标签位置在顶部
#     cbar.ax.get_yaxis().labelpad = 15
#     cbar.ax.set_position([0.83, 0.1, 0.05, 0.3])  # x, y, width, height
#     cbar.set_ticks([0,0.01, 0.02,0.03,0.04, 0.05])  # 设置固定的刻度值
#     cbar.ax.tick_params(labelsize=18)  # 设置刻度标签的字体大小

#     #cbar.ax.xaxis.set_ticks_position('top')  # 设置刻度也在顶部

#     # 添加基因数量图例
#     handles, labels = [], []
#     for size in [1, 2, 5]:  # 示例大小
#         handle = ax.scatter([], [], c='k', alpha=0.6, s=size* 300, label=f'{size}', linewidth=0)#edgecolors='w', 
#         handles.append(handle)
#     legend1 = ax.legend(
#     handles=handles,
#     title='Gene Count',
#     loc='upper left',
#     bbox_to_anchor=(0.95, 1.00),
#     fontsize='22',
#     title_fontsize='22',
#     frameon=False,
#     handlelength=1.5,  # 控制图例标记的长度
#     handletextpad=0.5,  # 增加标签文本和图例符号之间的间距
#     labelspacing=1  # 增加标签之间的垂直间距
# )

#     # 显示图形
#     plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/{title}.png', bbox_inches='tight')

# # 对三个数据集绘图
# process_and_plot(combined_df, 'Homo DEG GO function')


#==========================================================================================================================
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from icecream import install,ic
install()
config = {
    "font.family":'Times New Roman',
    "font.size": 25,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)
# 读取数据
speci = 'japonica'
bp_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/{speci}_BP_GO.csv')
cc_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/{speci}_CC_GO.csv')
mf_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/{speci}_MF_GO.csv')

combined_df = pd.concat([bp_data, cc_data, mf_data], ignore_index=True)

# 定义处理和绘图的函数
def process_and_plot(data, title):
    # 筛选Adjusted P-value小于0.05的行
    if speci == 'japonica':
        filtered_data = data[data['P-value'] < 0.05]
    else:
        filtered_data = data[data['Adjusted P-value'] < 0.05]
    filtered_data.to_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/{speci}_select_GO.csv',index=False)
    filtered_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_GO/{speci}_select_GO copy.csv')
    # 设置横坐标为Overlap中的分数，这里我们将它转换为整数分子以供后续使用
    richfactor = filtered_data['Overlap']
    richfactor_percent = richfactor.apply(lambda x: float(x.split('/')[0]) / float(x.split('/')[1]) * 100 if '/' in x else None)
    numerators = richfactor.apply(lambda x: int(x.split('/')[0]) if '/' in x and x.split('/')[0].isdigit() else None)

    # 设置基因数量表示气泡大小
    sizes = numerators

    # y轴设置为Term列的条目
    terms = filtered_data['Term']

    # 设置颜色映射Adjusted P-value
    if speci == 'japonica':
        colors = filtered_data['P-value']
    else:
        colors = filtered_data['Adjusted P-value']

    # 绘图
    if speci == 'homo':
        n=300
    elif speci =='indica' or speci =='zeamays':
        n=10
    elif speci =='japonica':
        n=50
    fig, ax = plt.subplots(figsize=(5, 7))
    scatter = ax.scatter(richfactor_percent, terms, s=sizes * n, c=colors, cmap='RdBu', edgecolors='black', linewidth=1, vmin=0, vmax=0.05,zorder=3) #alpha=0.6, 
    ax.set_title(title, fontsize=35, fontweight='bold',pad=10)
    ax.set_xlabel('Rich Factor (%)',fontsize=30,fontweight='bold') 
    ax.set_ylabel('GO Term',fontsize=30,fontweight='bold')
    ax.set_xticks([0, 25, 50, 75, 100])  # 设置x轴的固定刻度
    ax.set_xlim(-10, 110)
    #ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    # 绘制横向网格线
    ax.grid(True, axis='y', linestyle='--', linewidth=0.7, zorder=1)
    # 调整颜色条的位置，使其上下对半分

    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
    #cbar.set_label('FDR', rotation=0, labelpad=10)
    fig.text(0.77, 0.45, 'Q value', rotation=0, ha='left', va='center', fontsize=22)
    cbar.ax.xaxis.set_label_position('top')  # 设置标签位置在顶部
    cbar.ax.get_yaxis().labelpad = 15
    cbar.ax.set_position([0.83, 0.1, 0.05, 0.3])  # x, y, width, height
    cbar.set_ticks([0,0.01, 0.02,0.03,0.04, 0.05])  # 设置固定的刻度值
    cbar.ax.tick_params(labelsize=18)  # 设置刻度标签的字体大小

    #cbar.ax.xaxis.set_ticks_position('top')  # 设置刻度也在顶部
    if speci == 'homo':
        list = [1, 2, 5]
        labelspacing =1
    elif speci =='indica' or  speci =='zeamays':
        list = [5,10,30,60]
        labelspacing =0.6
    elif speci =='japonica':
        list = [5,10,30]
        labelspacing =0.6

    # 添加基因数量图例
    handles, labels = [], []
    for size in list:  # 示例大小
        handle = ax.scatter([], [], c='k', alpha=0.6, s=size* n, label=f'{size}', linewidth=0)#edgecolors='w', 
        handles.append(handle)
    legend1 = ax.legend(
    handles=handles,
    title='Gene Count',
    loc='upper left',
    bbox_to_anchor=(0.95, 1.00),
    fontsize='22',
    title_fontsize='22',
    frameon=False,
    handlelength=1.5,  # 控制图例标记的长度
    handletextpad=0.5,  # 增加标签文本和图例符号之间的间距
    labelspacing=labelspacing  # 增加标签之间的垂直间距
)

    # 显示图形
    plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/{title}.png', bbox_inches='tight')

# 对三个数据集绘图
process_and_plot(combined_df, f'{speci} DEG GO pathway')












