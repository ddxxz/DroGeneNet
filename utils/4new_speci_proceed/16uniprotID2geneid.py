
#通过uniprot的IDmapping网站将UniProtKB AC/ID转换成Ensembl Genomes
#==============================================================
# import pandas as pd

# # 读取TSV文件
# file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/idmapping_indicaBGI_initial.tsv'  # 替换为实际的文件路径
# data = pd.read_csv(file_path, sep='\t', header=0)

# # 筛选To列中以BGIO开头的行
# filtered_data = data[data['To'].str.startswith('BGIO')]

# # 保存筛选后的数据至CSV文件
# output_file_path = '/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/indica_BGI/idmapping_indicaBGI.tsv'  # 替换为实际的输出文件路径
# filtered_data.to_csv(output_file_path, index=False,sep='\t')

# print("Filtered data saved to:", output_file_path)
#=============================================================

import pandas as pd
for speci in ['wheat']:  #'indica_BGI','japonica','zeamays''sorghum'

    # 读取PPI数据
    ppi_file_path = f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/{speci}/PPI_reg.csv'  # 替换为实际的文件路径
    ppi_data = pd.read_csv(ppi_file_path)

    # 读取ID映射数据
    mapping_file_path = f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/{speci}/idmapping_{speci}.tsv'  # 替换为实际的文件路径
    mapping_data = pd.read_csv(mapping_file_path, sep='\t', header=None, names=['From', 'To'])

    protein_to_gene = dict(zip(mapping_data['From'], mapping_data['To']))

    # 替换PPI数据中的protein1和protein2列
    ppi_data['protein1'] = ppi_data['protein1'].map(protein_to_gene)
    ppi_data['protein2'] = ppi_data['protein2'].map(protein_to_gene)
    ppi_data = ppi_data.dropna()
    print(ppi_data)
    if speci=='sorghum':
        mapping_data1 = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.txt')
        genemap = dict(zip(mapping_data1['Gene stable ID'], mapping_data1['STRING ID']))
        ppi_data['protein1'] = ppi_data['protein1'].map(genemap)
        ppi_data['protein2'] = ppi_data['protein2'].map(genemap)
        ppi_data = ppi_data.dropna()
        ppi_data['protein1'] = ppi_data['protein1'].astype(str).str.split('.').str[1]
        ppi_data['protein2'] = ppi_data['protein2'].astype(str).str.split('.').str[1]
        print(ppi_data)

    # 删除含有NaN值的行（如果有映射失败的情况）
    ppi_data = ppi_data.dropna()
    print(ppi_data)

    # 保存转换后的数据至CSV文件
    output_file_path = f'/home/win/4T/GeneDL//data_output/GRN_result/Batch_correction/data/batch_correction_data/PPI/{speci}/{speci}_PPI_gene.pkl'  # 替换为实际的输出文件路径
    ppi_data.to_pickle(output_file_path)
    #ppi_data.to_csv(output_file_path, index=False)











