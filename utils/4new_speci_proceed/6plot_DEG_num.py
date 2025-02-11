
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from icecream import install,ic
install()
config = {
    "font.family":'Arial',
    "font.size": 18,
    "mathtext.fontset":'stix',
    "font.serif": ['Arial'],
}   #
plt.rcParams.update(config)

# 假设你的CSV文件名为"data.csv"
# file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/up_down_genenum.csv'
# data = pd.read_excel('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/DEG_homologues.xlsx',sheet_name='Sheet6')
# df=pd.DataFrame()
# # 
# deg_data_sheet2 = pd.read_excel('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/DEG_homologues.xlsx', sheet_name='Sheet6')
# deg_data_sheet1 = pd.read_excel('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_DEG/DEG_homologues.xlsx', sheet_name='Sheet5')

# # 获取Sheet2的列名列表
# columns = deg_data_sheet2.columns

# # 遍历每一列，找到None值并填充对应的Sheet1的值
# for col in columns:
#     # 检查当前列在Sheet1中是否存在
#     if col in deg_data_sheet1.columns:
#         # 进行填充：用Sheet1中相应列的值填充Sheet2中的NaN值
#         deg_data_sheet2[col] = deg_data_sheet2[col].fillna(deg_data_sheet1[col])
#     else:
#         print(f"Column {col} not found in Sheet1")
# # 使用pandas读取CSV文件
# #df = pd.read_csv(file_path,sep='\t')
# #colors = ['red', 'blue'] * 3 
# # 数据准备：将列名作为x轴标签，GeneNum作为y轴数据，并对索引进行重排以适应横置柱状图
# df['Group']=['indica_up','indica_down','japonica_up','japonica_down','zeamays_up','zeamays_down']
# data=deg_data_sheet2
# df['GeneNum']=[len(data['indica_up'].unique())-1,len(data['indica_down'].unique())-1,len(data['japonica_up'].unique())-1,len(data['japonica_down'].unique())-1,len(data['zeamays_up'].unique())-1,len(data['zeamays_down'].unique())-1]
data = {
    'Group': ['indica-up', 'indica-down', 'japonica-up', 'japonica-down', 'wheat-up', 'wheat-down','zeamays-up', 'zeamays-down','sorghum-up', 'sorghum-down',],
    'GeneNum': [943, 1479, 1928, 322,678,3391, 1531, 1165,854,1954]   #homo 86 59
}

# 创建 DataFrame
df = pd.DataFrame(data)

# 将 'Group' 列设置为索引
#df.set_index('Group', inplace=True)

data_to_plot = df[['Group', 'GeneNum']].set_index('Group')
ic(data_to_plot)
unique_groups = data_to_plot.index.unique()
ic(unique_groups)
colors = [(77,187,213), (77,187,213),(230,75,53),(230,75,53),(0,160,135),(0,160,135),(255, 205, 86),(255, 205, 86),(153, 102, 255),(153, 102, 255)]
colors_normalized = [(r/255, g/255, b/255) for r, g, b in colors]

# 绘制横置柱状图
# plt.figure(figsize=(15, 7))  # 设置图形大小

# # 使用颜色列表绘制柱状图
# ax = data_to_plot.plot(kind='barh', color=colors)
fig, ax = plt.subplots(figsize=(8, 5))
categories = df['Group'].values.tolist()
ic(categories)
# 直接遍历数据并为每个条形指定颜色

for index, value in enumerate(data_to_plot['GeneNum']):
    ax.barh(index, value, color=colors_normalized[index])

# 设置y轴的刻度位置
ax.set_yticks(range(len(data_to_plot['GeneNum'])))

# 然后设置刻度标签
ax.set_yticklabels(categories, fontsize=20, weight='bold')
# 在每个柱子上方添加数值标签
for i, (value, patch) in enumerate(zip(data_to_plot['GeneNum'], ax.patches)):
    # 计算文本位置，这里假设直接放在柱子顶部，可以根据需要调整偏移量
    text_x = patch.get_width() + 5  # 柱子右侧偏移一点
    text_y = patch.get_y() + patch.get_height() / 2  # 柱子中心高度
    ax.text(text_x, text_y, f'{value}', va='center', ha='left', fontsize=20)  # 添加文本，垂直居中，水平左对齐

# 添加标题和坐标轴标签

plt.title('Degree Up and Down',fontsize=30,weight='bold',pad=15)
plt.xlabel('Gene Degree',fontsize=25,weight='bold')
plt.ylabel('Groups',fontsize=25,weight='bold')
max_val = max(df['GeneNum']) * 1.5  # 增加10%的缓冲区
plt.xlim(0, max_val)
plt.gca().spines['top'].set_linewidth(1.5)
plt.gca().spines['bottom'].set_linewidth(1.5)
plt.gca().spines['right'].set_linewidth(1.5)
plt.gca().spines['left'].set_linewidth(1.5)
for label in ax.get_yticklabels():
    label.set_fontweight('bold')
for label in ax.get_xticklabels():
    label.set_fontweight('bold')
# 显示图形
plt.tight_layout()  # 自动调整子图参数, 使之填充整个图像区域
plt.savefig('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/deg_num.png',bbox_inches='tight')
























