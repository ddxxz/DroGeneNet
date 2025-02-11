
#-------------------------------获取Geneid和GOid---------------------------------------------------OSoutput
# import pandas as pd

# # 假设文件是以制表符分隔的，将 'sep' 参数替换成 ',' 如果是逗号分隔的CSV
# file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/Sbi_GO_annotation' # 或者 'Osj_GO_annotation.csv'
# output_file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/OSGOoutput.csv'

# # 读取文件，只获取前两列
# # header参数假设文件有标题行，如果没有请设定为None
# df = pd.read_csv(file_path, sep='\t', usecols=[0, 1,2], header=0)

# # 保存到新的CSV文件
# df.to_csv(output_file_path, index=False)

# print(f'文件已保存到 {output_file_path}')

#--------------------------------LOC转OS--------------------------------------------------------
import pandas as pd
from icecream import install,ic
install()
#===================================================
# 读取转换表，假设该表有两列: 'RAP' 和 'MSU'
# 提示: 根据转换表的实际情况修改列名
# conversion_table_path = '/mnt/g/GeneDL/OSDrought_GCNAT_Link/data_proceed/MSU2RAP/RAP-MSU_2018-03-29.txt'
# conversion_table = pd.read_csv(conversion_table_path)

# # 将转换表转换为字典，方便查询
# rap_to_msu_dict = pd.Series(conversion_table.iloc[:,1].values, index=conversion_table.iloc[:,0]).to_dict()

# relation={}
# for i in open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data_proceed/MSU2RAP/RAP-MSU_2018-03-29.txt"):
#     rap=str(i.split()[0])
#     msu=str(i.split()[1])
#     if msu!="None":
#         if "," in msu:
#             for a in msu.split(","):
#                 relation[a[0:-2]] = rap
#         else:
#             relation[msu[0:-2]] = rap
# ic(relation)
#===================================================

#===================================================
# relation = []
# # 使用open函数和with语句读取文件
# with open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt", 'r') as file:
#     for line in file:
#         # 分割每行的数据，只取前两个值（原始基因名和第一个映射基因名）
#         parts = line.strip().split('\t')
#         if len(parts) >= 2:  # 确保有足够的部分读取
#             original_gene = parts[0]
#             mapped_gene = parts[1]
#             # 保存结果
#             relation.append([original_gene, mapped_gene])

# # 将列表转换为DataFrame
# df3 = pd.DataFrame(relation, columns=['Original Gene', 'Mapped Gene'])
# # 创建映射字典
# gene_id_map = pd.Series(df3['Mapped Gene'].values, index=df3['Original Gene'].values).to_dict()
# # 假设你有一个DataFrame 'df' 包含需要转换的RAP标识符
# # 提示: 根据你的数据的实际情况调整代码
# #ic(gene_id_map)
# df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica_BGI/OSGOoutput.csv')

# #应用转换函数到DataFrame的列上
# df['OSGeneid'] = df['Gene_id'].map(gene_id_map)#df['Gene_id'] 

# #保存修改后的DataFrame
# df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica_BGI/OSGO.csv', index=False)
#===================================================
# map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_geneid_trans.csv')
# gene_id_map = dict(zip(map_txt['IWGSC1 CSS'], map_txt['IWGSC RefSeq v1.1']))
# df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/OSGOoutput.csv')
# df['OSGeneid'] = df['Gene_id'].map(gene_id_map)#df['Gene_id'] 
# df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/OSGO.csv', index=False)

#===================================================
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv')
gene_id_map = dict(zip(map_txt['merged_value'], map_txt['new_column']))
df = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/OSGOoutput.csv')
df['OSGeneid'] = df['Gene_id'].map(gene_id_map)#df['Gene_id'] 
df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/OSGO.csv', index=False)


#------------------------获取与基因表达矩阵一样的gene排列顺序以及对应的GOid---------------------------------------
# import pandas as pd
# from icecream import install,ic
# install()
# # 读取 all_genefile.csv 文件，它包含了 'Gene' 和 'Index' 列
# all_gene_file_path = '/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/ALL_Genefile.csv'
# all_gene_df = pd.read_csv(all_gene_file_path)

# # 读取 OSGO.csv 文件，它包含了 'OSGeneid' 和 'GO_ID' 列
# osgo_file_path = '/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/GOcluster/OSGO.csv'
# osgo_df = pd.read_csv(osgo_file_path)

# # 在all_gene_df中创建一个新列 'GO_ID' 并且默认填充为 'None'
# all_gene_df['GO_ID'] = 'None'

# # 创建 OSGeneid 到 GO_ID 的字典映射
# go_dict = pd.Series(osgo_df.GO_ID.values, index=osgo_df.OSGeneid).to_dict()

# # 对于 all_genefile.csv 中每行的 'Gene' 值，查找在 OSGO.csv 文件中的对应 'GO_ID'
# for index, row in all_gene_df.iterrows():
#     gene_id = row['Gene']
#     # 如果gene_id存在于go_dict字典中，则从字典中取得对应的GO_ID
#     if gene_id in go_dict:
#         all_gene_df.at[index, 'GO_ID'] = go_dict[gene_id]
#     # 如果gene_id不存在与go_dict字典中，默认已经是 'None'
# all_gene_df = all_gene_df.iloc[:, 1:]
# all_gene_df = pd.concat([all_gene_df.iloc[:, :1], all_gene_df.iloc[:, 2:]], axis=1)
# # 将更新后的 DataFrame 保存到新的csv文件中
# all_gene_df.to_csv('/mnt/g/GeneDL/OSDrought_GCNAT_Link/data/GOcluster/all_genefile_GO.csv', index=False)












