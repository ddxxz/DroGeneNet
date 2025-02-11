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
    "font.family":'Arial',
    "font.size": 22,
    "mathtext.fontset":'stix',
    "font.serif": ['Arial'],
    "font.weight": 'bold'
}
plt.rcParams.update(config)
# 替换'path_to_file.csv'为你的CSV文件的路径
titles = ['ck_cd','ck_ce']
for i,name in enumerate(['ck_cd','ck_ce']):

    indica_DEG_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/TF_indica_BGI_{name}_degree_exp_intersection.csv')
    japonica_DEG_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/TF_japonica_{name}_degree_exp_intersection.csv')
    zeamays_DEG_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/TF_zeamays_{name}_degree_exp_intersection.csv')
    wheat_DEG_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/TF_wheat_{name}_degree_exp_intersection.csv')
    sorghum_DEG_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_Degree_intersection_DEG/TF_sorghum_{name}_degree_exp_intersection.csv')


    indica_DEG = indica_DEG_data['0']
    japonica_DEG = japonica_DEG_data['0']
    zeamays_DEG = zeamays_DEG_data['0']
    wheat_DEG = wheat_DEG_data['0']
    sorghum_DEG = sorghum_DEG_data['0']

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

    map_txt2=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/homo/zeamay_geneid_trans.csv',sep='\t')
    gene_map2 = dict(zip(map_txt2['Gene stable ID'], map_txt2['Oryza sativa Japonica Group gene stable ID']))
    #ic(gene_map2)
    zeamays_DEG_map = zeamays_DEG.map(gene_map2)
    zeamays_DEG_map = zeamays_DEG_map.apply(lambda x: str(uuid.uuid4()) if pd.isna(x) else x)
    #ic(zeamays_DEG_map)

    map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_mapid.txt')
    gene_map = dict(zip(map_txt['Gene stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
    wheat_DEG_map = wheat_DEG.map(gene_map)
    wheat_DEG_map = wheat_DEG_map.apply(lambda x: str(uuid.uuid4()) if pd.isna(x) else x)

    map_txt=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_mapid.txt')
    gene_map = dict(zip(map_txt['Transcript stable ID'], map_txt['Oryza sativa Japonica Group gene stable ID']))
    sorghum_DEG_map = sorghum_DEG.map(gene_map)
    sorghum_DEG_map = sorghum_DEG_map.apply(lambda x: str(uuid.uuid4()) if pd.isna(x) else x)

    # japonica_DEG.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/japonica_DEG_{name}.csv',index=False)
    # indica_DEG_map.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/indica_DEG_map_{name}.csv',index=False)
    # zeamays_DEG_map.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/zeamays_DEG_map_{name}.csv',index=False)
        
    set1=set(japonica_DEG)
    set2=set(indica_DEG_map)
    set3=set(wheat_DEG_map)
    set4=set(zeamays_DEG_map)
    set5=set(sorghum_DEG_map)
    ic(set2)#set2,set3,set4,set5

    #colors = ['#8B0000', '#006400', '#00008B']

    labels = venn.get_labels([set1, set2, set3,set4,set5], fill=['number',
                                                        #'logic',,set3,set4
                                                    #'percent'
                                                            ]
                            )
    fig, ax = venn.venn5(labels,figsize=(6,6),fontsize=22, names=[f'japonica',f'indica',f'wheat',f'zeamays',f'sorghum',])#,f'leaf {speci} up',f'leaf {speci} down' f'root {speci} up',f'root {speci} down',dpi=96

    plt.title(f'homo {titles[i]}',fontsize=40,fontweight='bold')
    plt.savefig(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/paper_plot/new_framework_result/EXP_Degree_{name}.png',bbox_inches='tight')
    plt.clf()
    if name == 'DEG':
        common_gene= set1 & set2 & set3 & set4 & set5
        
        pd.DataFrame(common_gene).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/homo_deg.csv',index=False)

        japonica_individual = japonica_DEG
        indica_individual = indica_DEG
        wheat_individual = wheat_DEG
        zeamays_individual = zeamays_DEG
        sorghum_individual = sorghum_DEG

        pd.DataFrame(japonica_individual).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/japonica_deg.csv',index=False)
        pd.DataFrame(indica_individual).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/indica_deg.csv',index=False)
        pd.DataFrame(wheat_individual).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/wheat_deg.csv',index=False)
        pd.DataFrame(zeamays_individual).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/zeamays_deg.csv',index=False)
        pd.DataFrame(sorghum_individual).to_csv('/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_DEG/sorghum_deg.csv',index=False)






















