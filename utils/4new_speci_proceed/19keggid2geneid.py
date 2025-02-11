# import pandas as pd
# import os
# import glob
# from icecream import install,ic
# install()
# #=================================================获取转换数据==============================================================
# speci='wheat'
# # 定义文件夹路径
# folder_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/taes_edge'  # 替换为实际的文件夹路径

# # 获取所有txt文件路径
# txt_files = glob.glob(os.path.join(folder_path, '*.txt'))

# # 初始化一个空的DataFrame用于合并
# combined_df = pd.DataFrame()

# # 读取所有txt文件并合并
# for file in txt_files:
#     df = pd.read_csv(file, sep='\t', header=None, names=['protein1', 'protein2'])
#     combined_df = pd.concat([combined_df, df])

# # 保存合并后的数据至CSV文件
# output_combined_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/kegg_reg.csv'  # 替换为实际的输出文件路径
# combined_df.to_csv(output_combined_file_path, index=False)
# # ic(combined_df)
# # df_reversed = combined_df.rename(columns={"protein1": "protein2", "protein2": "protein1"})
# # # 将原始DataFrame与修改后的DataFrame进行合并，查找两个方向都存在的对
# # mutual_regulation = pd.merge(combined_df, df_reversed, on=["protein1", "protein2"])
# # # 查看结果
# # mutual_regulation_count = len(mutual_regulation) // 2  # 因为结果中每对会出现两次，所以需要除以2
# # ic(mutual_regulation_count)
# # 将两列拼成一列并获取唯一值
# unique_proteins = pd.unique(combined_df.values.ravel('K'))

# # 保存唯一蛋白质ID至另一个CSV文件
# unique_output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/kegg_unique_proteins.csv'  # 替换为实际的输出文件路径
# pd.DataFrame(unique_proteins, columns=['protein']).to_csv(unique_output_file_path, index=False)

# print("Combined data saved to:", output_combined_file_path)
# print("Unique proteins saved to:", unique_output_file_path)


#============================在uniprot网站，从KEGG转为UniprotKB，再从UniProtKB AC/ID转换成Ensembl Genomes==========================
'''
除小麦之外，其他几个物种的KEGG基因ID转换都可以在uniprot网站，从KEGG转为UniprotKB，再从UniProtKB AC/ID转换成Ensembl Genomes
小麦的NCBI Gene ID无法直接从uniprot网站上转，尝试多个网站均不行，突破口：从NCBI网站中找到GeneID对应的蛋白质ID，UniProtKB/TrEMBL，
再利用UniprotID转换为Ensembl Genomes
'''

# import pandas as pd

# # 读取TSV文件
# file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/idmapping_kegg_zeamays.tsv'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep='\t', header=None)

# # 筛选前两列
# data_filtered = data.iloc[:, :2]

# # 保存筛选后的数据至CSV文件
# output_file_path = '/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zeamays/kegg_pathway_proteinid.csv'  # 替换为实际的输出文件路径
# data_filtered.to_csv(output_file_path, index=False, header=False)

# print("Filtered data saved to:", output_file_path)

#===============================================================================

import pandas as pd
from icecream import install,ic
install()
for speci in ['wheat']: #'sorghum''japonica','zeamays'

    # 读取PPI数据
    ppi_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/kegg_reg.csv'  # 替换为实际的文件路径
    ppi_data = pd.read_csv(ppi_file_path)
    ic(ppi_data)
    # 读取ID映射数据

    if speci=='sorghum':
        mapping_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/sorghum/keggid_geneid.csv'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path, names=['From', 'Gene Names'],sep='\t')
        protein_to_gene1 = dict(zip(mapping_data['From'], mapping_data['Gene Names']))

        mapping_file_path = f'/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.txt'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path)
        protein_to_gene2 = dict(zip(mapping_data['Gene stable ID'], mapping_data['STRING ID']))
    elif speci=='wheat':
        mapping_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/geneid2uniprot.csv'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path, names=['GeneID', 'UniProtKB/TrEMBL'])
        protein_to_gene1 = dict(zip(mapping_data['GeneID'], mapping_data['UniProtKB/TrEMBL']))
        #ic(protein_to_gene1)
        mapping_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/idmapping_wheat.tsv'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path,sep='\t')
        protein_to_gene2 = dict(zip(mapping_data['From'], mapping_data['To']))
    else:
        mapping_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/kegg_pathway_proteinid.csv'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path, header=None, names=['From', 'Entry'])
        protein_to_gene1 = dict(zip(mapping_data['From'], mapping_data['Entry']))
        ic(protein_to_gene1)

        mapping_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/idmapping_kegg_{speci}.tsv'  # 替换为实际的文件路径
        mapping_data = pd.read_csv(mapping_file_path, sep='\t', header=None, names=['From', 'To'])
        protein_to_gene2 = dict(zip(mapping_data['From'], mapping_data['To']))

    # 替换PPI数据中的protein1和protein2列
    ic(ppi_data)
    if speci =='wheat':
        #ic(ppi_data['protein1'].astype(str).str.split(':').str[1])
        ppi_data['protein1'] = ppi_data['protein1'].astype(str).str.split(':').str[1].map(protein_to_gene1)
        ppi_data['protein2'] = ppi_data['protein2'].astype(str).str.split(':').str[1].map(protein_to_gene1)
    else:
        ppi_data['protein1'] = ppi_data['protein1'].map(protein_to_gene1)
        ppi_data['protein2'] = ppi_data['protein2'].map(protein_to_gene1)
    ic(ppi_data)
    ppi_data['protein1'] = ppi_data['protein1'].map(protein_to_gene2)
    ppi_data['protein2'] = ppi_data['protein2'].map(protein_to_gene2)
    # 删除含有NaN值的行（如果有映射失败的情况）
    #ic(protein_to_gene2)
    ic(ppi_data)
    ppi_data = ppi_data.dropna()
    if speci=='sorghum':
        ppi_data['protein1'] = ppi_data['protein1'].astype(str).str.split('.').str[1]
        ppi_data['protein2'] = ppi_data['protein2'].astype(str).str.split('.').str[1]
    ic(ppi_data)
    # 保存转换后的数据至CSV文件
    output_file_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/{speci}/{speci}_kegg_gene.csv'  # 替换为实际的输出文件路径
    #ppi_data.to_pickle(output_file_path)
    ppi_data.to_csv(output_file_path, index=False)

    print("Mapped PPI data saved to:", output_file_path)



