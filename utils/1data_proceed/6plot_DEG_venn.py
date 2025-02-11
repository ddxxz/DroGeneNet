import pandas as pd
from icecream import install,ic
install()
import matplotlib.pyplot as plt
import venn
from matplotlib_venn import venn3
import pandas as pd
import string
import uuid
config = {
    "font.family":'Times New Roman',
    "font.size": 22,
    "mathtext.fontset":'stix',
    "font.serif": ['Times New Roman'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)
# 替换'path_to_file.csv'为你的CSV文件的路径
titles = ['DEG','DEG up','DEG down']
for i,name in enumerate(['DEG','DEG_up','DEG_down']):

    indica_DEG_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/CK_CD/indica_BGI/indica_BGI_CK_CD_{name}.csv')
    japonica_DEG_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/CK_CD/japonica/japonica_CK_CD_{name}.csv')
    zeamays_DEG_data = pd.read_csv(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/CK_CD/zeamays/zeamays_CK_CD_{name}.csv')

    indica_DEG = indica_DEG_data['FEATURE_ID']
    japonica_DEG = japonica_DEG_data['FEATURE_ID']
    zeamays_DEG = zeamays_DEG_data['FEATURE_ID']
    #ic(zeamays_DEG)

    map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/indica_gene_mapos.txt')
    gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
    indica_DEG_map = indica_DEG.map(gene_map)
    indica_DEG_map = indica_DEG_map.apply(lambda x: str(uuid.uuid4()) if pd.isna(x) else x)
    #ic(indica_DEG_map)

    file_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt'
    map_txt1 = pd.read_csv(
        file_path,
        sep='\t',
        usecols=[0, 1],  # 只读取前两列
        names=['input', '0'],  # 假设文件没有表头，手动指定列名
        header=None  # 告诉pandas没有表头
    )
    #map_txt1=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt',sep='\t')
    # gene_map1 = dict(zip(map_txt1['input'], map_txt1['0']))
    # zeamays_DEG1 = zeamays_DEG.map(gene_map1)
    #ic(gene_map1)
    map_txt2=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')
    gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
    #ic(gene_map2)
    zeamays_DEG_map = zeamays_DEG.map(gene_map2)
    zeamays_DEG_map = zeamays_DEG_map.apply(lambda x: str(uuid.uuid4()) if pd.isna(x) else x)
    #ic(zeamays_DEG_map)

    # japonica_DEG.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/japonica_DEG_{name}.csv',index=False)
    # indica_DEG_map.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/indica_DEG_map_{name}.csv',index=False)
    # zeamays_DEG_map.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/zeamays_DEG_map_{name}.csv',index=False)
        
    set1=set(japonica_DEG)
    set2=set(indica_DEG_map)
    set3=set(zeamays_DEG_map)

    colors = ['#8B0000', '#006400', '#00008B']

    labels = venn.get_labels([set1, set2, set3], fill=['number',
                                                        #'logic',,set3,set4
                                                    #'percent'
                                                            ]
                            )
    fig, ax = venn.venn3(labels,figsize=(6,6),fontsize=22, names=[f'japonica',f'indica',f'zeamays'])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

    plt.title(f'homo {titles[i]}',fontsize=40,fontweight='bold')
    plt.savefig(f'/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/homo_{name}.png',bbox_inches='tight')
    plt.clf()
    if name == 'DEG':
        common_gene= set1 & set2 & set3
        
        pd.DataFrame(common_gene).to_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/homo_deg.csv',index=False)

        japonica_individual = japonica_DEG
        indica_individual = indica_DEG
        zeamays_individual = zeamays_DEG
        pd.DataFrame(japonica_individual).to_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/japonica_deg.csv',index=False)
        pd.DataFrame(indica_individual).to_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/indica_deg.csv',index=False)
        pd.DataFrame(zeamays_individual).to_csv('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/zeamays_deg.csv',index=False)






















