# import subprocess
# import pandas as pd
# from icecream import install,ic
# install()
# from bs4 import BeautifulSoup

# data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/kegg_unique_proteins.csv')
# genes = data['protein'][:10].str.split(':').str[1]


# import requests

# # 请求结果字典
# results = {}

# # 循环查询每个基因
# for gene in genes:
#     url = f"https://www.ricedata.cn/gene/accessions_switch.aspx?para={gene}&genenm=&cloned=false&located=false&chro="
#     response = requests.get(url)  # 发送GET请求
    
#     # 确认请求成功
#     if response.status_code == 200:
#         # 这里假设返回数据是HTML，你也许需要使用 Beautiful Soup 处理它
#         results[gene] = response.text  # 或者 response.json() 如果返回是JSON
#     else:
#         results[gene] = "Error or no data"

# # 打印或处理结果
# all_data =[]
# for gene, result in results.items():
#     #print(f"Gene: {gene}: Result: {result}")

#     # 如果需要解析HTML，你可能要安装并导入BeautifulSoup来解析response.text
    
#     # HTML 数据
#     html_data = result# """
#     #您的HTML串在这里...
#     #"""

#     soup = BeautifulSoup(html_data, 'html.parser')

#     # 查找所有<tr>标签，每个<tr>标签代表一行
#     rows = soup.find_all('tr', bgcolor='#FFFFFF')

#     # 用于保存解析后的数据
#     data = []

#     for row in rows:
#         # 在当前行中查找所有<td>标签，每个<td>标签是一个单元格
#         cells = row.find_all('td', style='border-bottom:1px solid silver')
        
#         # 从每个单元格中提取文本(如果需要的话，还可以提取<a>标签的href属性等)
#         row_data = [cell.get_text(strip=True) for cell in cells]
        
#         # 将提取的行数据添加到数据列表中
#         data.append(row_data)
#     all_data.append(data[2])

# ic(all_data)

# all_data_df = pd.DataFrame(all_data)
# all_data_df.columns = ['GeneID','基因名称或注释','基因符号','RAP_Locus','MSU_Locus或其它','NCBI_Locus','cDNAs','RefSeq_Locus','Uniprots']
# ic(all_data_df)
# all_data_df.to_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/outdata/CK_CE_DEG_gene_function_name.csv',index=False)

import requests
import pandas as pd
from bs4 import BeautifulSoup


def fetch_uniprot_id(geneid):
    # 构建查询URL
    url = f"https://www.ncbi.nlm.nih.gov/gene/?term={geneid}"

    # 发出请求
    response = requests.get(url)
    
    # 检查请求是否成功
    if response.status_code != 200:
        print(f"Failed to retrieve data for GeneID {geneid}")
        return None

    # 使用BeautifulSoup解析HTML内容
    soup = BeautifulSoup(response.text, 'html.parser')

    # 查找包含UniProtKB/TrEMBL的dt标签
    uniprot_dt = soup.find('dt', string='UniProtKB/TrEMBL')
    #print('UniProtKB/TrEMBL')
    if uniprot_dt:
        # 找到相邻的dd标签中的链接
        uniprot_dd = uniprot_dt.find_next_sibling('dd')
        if uniprot_dd:
            uniprot_link = uniprot_dd.find('a')
            if uniprot_link:
                # 提取UniProt ID
                uniprot_id = uniprot_link.text.strip()
                return uniprot_id
    else:
        print(f"UniProtKB/TrEMBL section not found for GeneID {geneid}")
        return None

# 测试函数，使用指定的GeneID
# geneid = "123096499"
# uniprot_id = fetch_uniprot_id(geneid)
# if uniprot_id:
#     print(f"UniProtKB/TrEMBL ID for GeneID {geneid}: {uniprot_id}")

data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/kegg_unique_proteins.csv')
genes = data['protein'].str.split(':').str[1]

result_df = pd.DataFrame(columns=['GeneID', 'UniProtKB/TrEMBL'])
# 遍历基因列表，获取每个基因的 UniProt ID
for geneid in genes:
    uniprot_id = fetch_uniprot_id(geneid)
    # 将结果保存到 DataFrame 中
    result_df = result_df.append({'GeneID': geneid, 'UniProtKB/TrEMBL': uniprot_id}, ignore_index=True)

# 打印或保存结果
print(result_df)
result_df.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/geneid2uniprot.csv',index=False)




