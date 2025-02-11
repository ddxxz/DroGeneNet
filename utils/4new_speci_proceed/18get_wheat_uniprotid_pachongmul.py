import requests
import pandas as pd
from bs4 import BeautifulSoup
from concurrent.futures import ThreadPoolExecutor, as_completed
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from tqdm import tqdm

# 创建全局会话
session = requests.Session()

# 配置重试机制，提升请求稳定性
retries = Retry(total=5, backoff_factor=0.1, status_forcelist=[500, 502, 503, 504])
session.mount('https://', HTTPAdapter(max_retries=retries))

# 定义函数获取 UniProt ID
def fetch_uniprot_id(geneid):
    # 构建查询URL
    url = f"https://www.ncbi.nlm.nih.gov/gene/?term={geneid}"

    try:
        # 发出请求
        response = session.get(url, timeout=10)
        
        # 检查请求是否成功
        if response.status_code != 200:
            print(f"Failed to retrieve data for GeneID {geneid}")
            return geneid, None

        # 使用BeautifulSoup解析HTML内容
        soup = BeautifulSoup(response.text, 'html.parser')

        # 查找包含UniProtKB/TrEMBL的dt标签
        uniprot_dt = soup.find('dt', string='UniProtKB/TrEMBL')

        if uniprot_dt:
            # 找到相邻的dd标签中的链接
            uniprot_dd = uniprot_dt.find_next_sibling('dd')
            if uniprot_dd:
                uniprot_link = uniprot_dd.find('a')
                if uniprot_link:
                    # 提取UniProt ID
                    uniprot_id = uniprot_link.text.strip()
                    return geneid, uniprot_id

        return geneid, None

    except Exception as e:
        print(f"Error fetching data for GeneID {geneid}: {e}")
        return geneid, None

# 读取CSV文件
data = pd.read_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/kegg_unique_proteins.csv')
genes = data['protein'].str.split(':').str[1]

# 线程数，可以根据你的机器性能调整。过高可能会导致请求失败
MAX_WORKERS = 30

# 并行处理函数，带进度条
def fetch_uniprot_ids_parallel(genes):
    results = []
    
    # 使用 ThreadPoolExecutor 进行并行处理
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        # 提交任务并发执行
        futures = [executor.submit(fetch_uniprot_id, geneid) for geneid in genes]
        
        # 使用 tqdm 包裹进度条
        for future in tqdm(as_completed(futures), total=len(futures), desc="Fetching UniProt IDs"):
            geneid, uniprot_id = future.result()
            results.append({'GeneID': geneid, 'UniProtKB/TrEMBL': uniprot_id})
    
    return results

# 进行并行爬取
results = fetch_uniprot_ids_parallel(genes)

# 将结果转换为 DataFrame
result_df = pd.DataFrame(results)

# 打印或保存结果
print(result_df)
result_df.to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/wheat/geneid2uniprotnew.csv', index=False)
