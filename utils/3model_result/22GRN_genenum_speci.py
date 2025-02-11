# 自己构建的调控网络的基因数量

import pandas as pd
from icecream import install,ic
install()
import pickle

# 标签的调控网络的基因数量
#=================================================indica  label  count===================================================
file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_regulation_merged_Osi.txt'
data = pd.read_csv(file_path, sep='\t')

# 合并 TF 列和 Target 列
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,2]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_gene_mapos.txt')
indica_gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
gene_counts['OSgene'] = gene_counts['Gene'].map(indica_gene_map)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/TF_indica_BGI_label_gene_counts.csv', index=False)


with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/indica_BGI_PPI_gene.pkl','rb') as fp:
    data = pickle.load(fp)
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,1]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_gene_mapos.txt')
indica_gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
gene_counts['OSgene'] = gene_counts['Gene'].map(indica_gene_map)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/PPI_indica_BGI_label_gene_counts.csv', index=False)

data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica/indica_kegg_gene.csv')

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/indica_KEGG_gene.pkl','rb') as fp:
#     data = pickle.load(fp)
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,1]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_gene_mapos.txt')
indica_gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
gene_counts['OSgene'] = gene_counts['Gene'].map(indica_gene_map)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/KEGG_indica_label_gene_counts.csv', index=False)


#=================================================japonica  label  count===================================================
file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/japonica/japonica_regulation_merged_Osj.txt'
data = pd.read_csv(file_path, sep='\t')
# 合并 TF 列和 Target 列
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,2]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
relation={}
for i in open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data_proceed/MSU2RAP/RAP-MSU_2018-03-29.txt"):
    rap=str(i.split()[0])
    msu=str(i.split()[1])
    if msu!="None":
        if "," in msu:
            for a in msu.split(","):
                relation[a[0:-2]] = rap
        else:
            relation[msu[0:-2]] = rap
def convert_rap_to_msu(rap_id):
    return relation.get(rap_id, rap_id)
def convert_msu_to_rap(msu_id):
    return relation.get(msu_id, msu_id)
gene_counts['OSgene'] = gene_counts['Gene'].apply(convert_msu_to_rap)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/TF_japonica_label_gene_counts.csv', index=False)

with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica_PPI_gene.pkl','rb') as fp:
    data = pickle.load(fp)
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,1]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
gene_counts['OSgene'] = gene_counts['Gene']#.apply(convert_msu_to_rap)
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/PPI_japonica_label_gene_counts.csv', index=False)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/japonica/japonica_KEGG_gene.pkl','rb') as fp:
#     data = pickle.load(fp)
data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/japonica_kegg_gene.csv')
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,1]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
gene_counts['OSgene'] = gene_counts['Gene']#.apply(convert_msu_to_rap)
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/KEGG_japonica_label_gene_counts.csv', index=False)


#=================================================zeamays label  count===================================================
file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/regulation_merged_Zma.txt'
data = pd.read_csv(file_path, sep='\t')
# 合并 TF 列和 Target 列
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,2]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt'
map_txt1 = pd.read_csv(
    file_path,
    sep='\t',
    usecols=[0, 1],  # 只读取前两列
    names=['input', '0'],  # 假设文件没有表头，手动指定列名
    header=None  # 告诉pandas没有表头
)
#map_txt1=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt',sep='\t')
zm_gene_map1 = dict(zip(map_txt1['input'], map_txt1['0']))
gene_counts['ZMgene'] = gene_counts['Gene'].map(zm_gene_map1)
map_txt2=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')
zm_gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
gene_counts['OSgene'] = gene_counts['ZMgene'].map(zm_gene_map2)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/TF_zeamays_label_gene_counts.csv', index=False)


with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/zeamays/zeamays_PPI_gene.pkl','rb') as fp:
    data = pickle.load(fp)

merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,2]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
gene_counts['OSgene'] = gene_counts['Gene'].map(zm_gene_map2)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/PPI_zeamays_label_gene_counts.csv', index=False)

# with open('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_KEGG_gene.pkl','rb') as fp:
#     data = pickle.load(fp)
data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/zeamays_kegg_gene.csv')
merged_genes = pd.concat([data.iloc[:,0], data.iloc[:,1]])
# 统计每个基因名出现的次数
gene_counts = merged_genes.value_counts().reset_index()
gene_counts.columns = ['Gene', 'Count']
gene_counts['OSgene'] = gene_counts['Gene'].map(zm_gene_map2)
# 显示结果
ic(gene_counts)
# 如果需要保存结果到新的 CSV 文件
gene_counts.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/KEGG_zeamays_label_gene_counts.csv', index=False)


#================================水分条件变化下的基因数量==========================================
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type== 'KEGG':
        species=['indica','zeamays','japonica']
    else:
        species=['indica_BGI','zeamays','japonica']
    for speci in species:
        for treat in ['CK','CD','CE']:
            # 读取 txt 文件，假设文件是制表符分隔的
            file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/model_pred_edgeweight/{data_type}_{speci}/{treat}_all.txt'
            data = pd.read_csv(file_path, sep='\t')

            # 合并 TF 列和 Target 列
            merged_genes = pd.concat([
                data[['TF', 'TF_GO']].rename(columns={'TF': 'Gene', 'TF_GO': 'GO'}),
                data[['Target', 'Target_GO']].rename(columns={'Target': 'Gene', 'Target_GO': 'GO'})
            ])
            # 统计每个基因名出现的次数
            gene_counts = merged_genes.value_counts().reset_index()
            gene_counts.columns = ['Gene','GO', 'Count']
            if speci =='indica_BGI':
                gene_counts['OSgene'] = gene_counts['Gene'].map(indica_gene_map)
            elif speci =='zeamays':
                gene_counts['OSgene'] = gene_counts['Gene'].map(zm_gene_map2)
            else:
                gene_counts['OSgene'] = gene_counts['Gene']
            # 显示结果
            ic(gene_counts)

            # 如果需要保存结果到新的 CSV 文件
            gene_counts.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_{speci}_{treat}_gene_counts.csv', index=False)

#===============================================================================================

import pandas as pd
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    for treat in ['CK','CD','CE','label']:
        # 读取三个 CSV 文件
        if data_type=='KEGG':
            file_paths = [
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_indica_{treat}_gene_counts.csv',
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_japonica_{treat}_gene_counts.csv',
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_zeamays_{treat}_gene_counts.csv'
            ]   
        else:         
            file_paths = [
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_indica_BGI_{treat}_gene_counts.csv',
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_japonica_{treat}_gene_counts.csv',
                f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_zeamays_{treat}_gene_counts.csv'
            ]

        dfs = [pd.read_csv(file) for file in file_paths]

        # 重命名每个数据框中的 Count 列
        dfs[0] = dfs[0].rename(columns={'Count': 'indica_Count'})
        dfs[1] = dfs[1].rename(columns={'Count': 'japonica_Count'})
        dfs[2] = dfs[2].rename(columns={'Count': 'zeamays_Count'})
        
        for df in dfs:
            df['OSgene'] = df['OSgene'].astype(str)
        # 合并数据框，确保相同的 OSgene 行拼接在一起
        merged_df = dfs[0].merge(dfs[1], on='OSgene', how='inner').merge(dfs[2], on='OSgene', how='inner')

        final_df = merged_df[['OSgene', 'indica_Count', 'japonica_Count', 'zeamays_Count']]
        
        if final_df.duplicated(subset=['OSgene']).any():
            # 找到重复的 OSgene
            duplicates = final_df[final_df.duplicated(subset=['OSgene'], keep=False)]
            
            # 保留前两列的第一个值，并累加最后一列的 Count
            unique_df = duplicates.groupby('OSgene').agg({
                'indica_Count': 'first',
                'japonica_Count': 'first',
                'zeamays_Count': 'first'
            }).reset_index()
            
            # 去除重复的行
            final_df = final_df.drop_duplicates(subset=['OSgene'], keep=False)
            
            # 将累加后的行添加回数据框
            final_df = pd.concat([final_df, unique_df], ignore_index=True)

        print(merged_df)

        # 保存结果到新的 CSV 文件
        output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_{treat}_genenum.csv'
        final_df.to_csv(output_file_path, index=False)



#==============================得到靶基因调控网络差异大的行===================================

import pandas as pd

# 读取 CSV 文件
file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_label_genenum.csv'
data = pd.read_csv(file_path)

filtered_data = data[(data.iloc[:, 1:] < 100).all(axis=1)]

# 计算第一二列求平均后与第三列的差值
filtered_data['Diff'] = (filtered_data.iloc[:, 1] + filtered_data.iloc[:, 2]) / 2 - filtered_data.iloc[:, 3]

# 计算差值的绝对值
filtered_data['AbsDiff'] = filtered_data['Diff'].abs()

# 根据差值的绝对值排序
sorted_data = filtered_data.sort_values(by='AbsDiff',ascending=False)

# 显示结果
ic(sorted_data)

# 如果需要保存结果到新的 CSV 文件
output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_target_gene_labeldata.csv'
sorted_data.to_csv(output_file_path, index=False)





go_terms = ['GO:0009399', 'GO:0006808', 'GO:0071705', 'GO:0006807', 'GO:0019740','GO:0009414']

# 将 GO 术语列表转换为正则表达式
go_pattern = '|'.join(go_terms)

data_types=['TF','PPI','KEGG']
for data_type in data_types:
    for treat in ['CK','CD','CE','label']:
        file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_{treat}_genenum.csv'
        data = pd.read_csv(file_path)

        go_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_japonica_CK_gene_counts.csv')
        
        go_data = go_data[['OSgene', 'GO']]
        # 将go_data中的GO列合并到data中
        data = data.merge(go_data, on='OSgene', how='left')
        ic(data.shape)
        filtered_data = data#[(data.iloc[:, 1:] > 500).all(axis=1)]

        # 计算第一二列求平均后与第三列的差值
        filtered_data['Diff'] = (filtered_data.iloc[:, 1] + filtered_data.iloc[:, 2]) / 2 - filtered_data.iloc[:, 3]

        # 计算差值的绝对值
        filtered_data['AbsDiff'] = filtered_data['Diff'].abs()

        # 根据差值的绝对值排序
        sorted_data = filtered_data.sort_values(by='AbsDiff',ascending=False)

        sorted_data = sorted_data.dropna(subset=['GO'])
        # 筛选GO列中包含指定GO术语的行
        sorted_data = sorted_data[sorted_data['GO'].str.contains(go_pattern, regex=True)]

        # 显示结果
        ic(sorted_data)

        # 如果需要保存结果到新的 CSV 文件
        output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_sort_nitrogen_oryza_vs_zm_{treat}data.csv'
        sorted_data.to_csv(output_file_path, index=False)



import pandas as pd
from icecream import ic
data_types=['TF','PPI','KEGG']
for data_type in data_types:
    for treat in ['CK', 'CD', 'CE', 'label']:
        file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_{treat}_genenum.csv'
        data = pd.read_csv(file_path)
        
        go_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_japonica_CK_gene_counts.csv')
         
        go_data = go_data[['OSgene', 'GO']]
        # 将go_data中的GO列合并到data中
        data = data.merge(go_data, on='OSgene', how='left')

        filtered_data = data  # 如果有筛选条件，可以取消注释并修改

        # 计算第一列和第二列的差值
        filtered_data['Diff'] = filtered_data.iloc[:, 1] - filtered_data.iloc[:, 2]

        # 计算差值的绝对值
        filtered_data['AbsDiff'] = filtered_data['Diff'].abs()

        # 根据差值的绝对值排序
        sorted_data = filtered_data.sort_values(by='AbsDiff', ascending=False)
        sorted_data = sorted_data.dropna(subset=['GO'])
        # 筛选GO列中包含指定GO术语的行
        sorted_data = sorted_data[sorted_data['GO'].str.contains(go_pattern, regex=True)]
        # 显示结果
        ic(sorted_data)

        # 如果需要保存结果到新的 CSV 文件
        output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_genenum/{data_type}_sort_nitrogen_indica_vs_japonica_{treat}data.csv'
        sorted_data.to_csv(output_file_path, index=False)








