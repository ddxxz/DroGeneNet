import numpy as np
import pandas as pd

data_types=['TF','PPI','KEGG']
for data_type in data_types:
    if data_type== 'KEGG':
        species=['wheat','sorghum']#'indica','zeamays','japonica',
    else:
        species=['wheat','sorghum']#'indica_BGI','zeamays','japonica',
    for speci in species:
        if speci=='zeamays':
            relation={}
            for i in open("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/gene_id_trans.txt"):
                rap=str(i.split()[0])
                msu=str(i.split()[1])
                #ic(rap)
                #ic(msu)
                if msu!="None":
                    if " " in msu:
                        for a in msu.split(" "):
                            relation[a] = rap
                    else:
                        relation[msu] = rap
            #ic(relation)
            # 读取目标 CSV 文件中的某一列为 DataFrame
            input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/zeamays/Zma_TF_list.txt'  # 替换为真实的 csv 文件路径
            df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
            def invert_dict(original_dict):
                """翻转字典，使得原始的值变成键，原始的键变成值。"""
                inverted_dict = {}
                for key, value in original_dict.items():
                    if value in inverted_dict:
                        # 确保字典的值是唯一的，如果不唯一，则添加到列表中
                        if not isinstance(inverted_dict[value], list):
                            inverted_dict[value] = [inverted_dict[value]]
                        inverted_dict[value].append(key)
                    else:
                        inverted_dict[value] = key
                return inverted_dict

            def map_values_to_keys(value_list, mapping_dict):
                """使用翻转后的字典映射列表中的值到对应的键。"""
                inverted_dict = invert_dict(mapping_dict)
                return [inverted_dict.get(value, None) for value in value_list]
            data1 = np.array(df.iloc[:, 1].apply(lambda x: f'{x}'))
            #ic(data1)
            mapped_keys = map_values_to_keys(data1,relation)
            #ic(mapped_keys)
            #ic(df.iloc[:,1])
            # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
            #keys_found = [key for key, val in relation.items() if val in df.iloc[:,1]]
            df['mapid'] = mapped_keys#df.iloc[:,1].map(relation).fillna("None")
            #ic(mapped_keys)
            #ic(df['Os_geneid'])
            df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_mulgeneid/{data_type}_zm_tf_list.csv',index=False)
            ck_ce_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_zeamays/CK_CE/all_zeamays_CK_CE_DEG.csv')
            ck_ce_degree_gene = ck_ce_degree['Unnamed: 0'][:10]
            mask = ck_ce_degree_gene.isin(df['mapid'])
            existing_values = ck_ce_degree_gene[mask]
            print(existing_values)

            ck_cd_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_zeamays/CK_CD/all_zeamays_CK_CD_DEG.csv')
            ck_cd_degree_gene = ck_cd_degree['Unnamed: 0'][:10]
            mask = ck_cd_degree_gene.isin(df['mapid'])
            existing_values = ck_cd_degree_gene[mask]
            print(existing_values)

            ck_ce_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_zeamays/CK_CE/zeamays_CK_CE_DEG.csv')
            ck_ce_exp_gene = ck_ce_exp['FEATURE_ID'][:10]
            mask = ck_ce_exp_gene.isin(df['mapid'])
            existing_values = ck_ce_exp_gene[mask]
            print(existing_values)

            ck_cd_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_zeamays/CK_CD/zeamays_CK_CD_DEG.csv')
            ck_cd_exp_gene = ck_cd_exp['FEATURE_ID'][:10]
            mask = ck_cd_exp_gene.isin(df['mapid'])
            existing_values = ck_cd_exp_gene[mask]
            print(existing_values)

        elif speci=='indica_BGI':

            input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/Indica_Osi_TF_list.txt'  # 替换为真实的 csv 文件路径
            df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
            #ic(df)
            # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
            df['mapid'] = df.iloc[:,1]
            df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_mulgeneid/{data_type}_indica_tf_list.csv',index=False)
            
            ck_ce_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_indica_BGI/CK_CE/all_indica_BGI_CK_CE_DEG.csv')
            ck_ce_degree_gene = ck_ce_degree['Unnamed: 0'][:10]
            mask = ck_ce_degree_gene.isin(df['mapid'])
            existing_values = ck_ce_degree_gene[mask]
            print(existing_values)

            ck_cd_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_indica_BGI/CK_CD/all_indica_BGI_CK_CD_DEG.csv')
            ck_cd_degree_gene = ck_cd_degree['Unnamed: 0'][:10]
            mask = ck_cd_degree_gene.isin(df['mapid'])
            existing_values = ck_cd_degree_gene[mask]
            print(existing_values)

            ck_ce_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_indica_BGI/CK_CE/indica_BGI_CK_CE_DEG.csv')
            ck_ce_exp_gene = ck_ce_exp['FEATURE_ID'][:10]
            mask = ck_ce_exp_gene.isin(df['mapid'])
            existing_values = ck_ce_exp_gene[mask]
            print(existing_values)

            ck_cd_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_indica_BGI/CK_CD/indica_BGI_CK_CD_DEG.csv')
            ck_cd_exp_gene = ck_cd_exp['FEATURE_ID'][:10]
            mask = ck_cd_exp_gene.isin(df['mapid'])
            existing_values = ck_cd_exp_gene[mask]
            print(existing_values)
        
        elif speci=='wheat':

            input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/Tae_TF_list.txt'  # 替换为真实的 csv 文件路径
            df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
            mapdata = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/wheat/wheat_geneid_trans.csv')
            gene_id_map = pd.Series(mapdata['IWGSC RefSeq v1.1'].values, index=mapdata['IWGSC1 CSS'].values).to_dict()
            df['Gene_ID'] = df['Gene_ID'].map(gene_id_map)
            #ic(df)
            # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
            df['mapid'] = df.iloc[:,1]
            print(df['mapid'])
            df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_mulgeneid/{data_type}_wheat_tf_list.csv',index=False)
            
            ck_ce_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_wheat/CK_CE/all_wheat_CK_CE_DEG.csv')
            ck_ce_degree_gene = ck_ce_degree['Unnamed: 0']#[:10]
            mask = ck_ce_degree_gene.isin(df['mapid'])
            existing_values = ck_ce_degree_gene[mask]
            print(existing_values)

            ck_cd_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_wheat/CK_CD/all_wheat_CK_CD_DEG.csv')
            ck_cd_degree_gene = ck_cd_degree['Unnamed: 0']#[:10]
            mask = ck_cd_degree_gene.isin(df['mapid'])
            existing_values = ck_cd_degree_gene[mask]
            print(existing_values)

            ck_ce_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_wheat/CK_CE/wheat_CK_CE_DEG.csv')
            ck_ce_exp_gene = ck_ce_exp['FEATURE_ID']#[:10]
            mask = ck_ce_exp_gene.isin(df['mapid'])
            existing_values = ck_ce_exp_gene[mask]
            print(existing_values)

            ck_cd_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_wheat/CK_CD/wheat_CK_CD_DEG.csv')
            ck_cd_exp_gene = ck_cd_exp['FEATURE_ID']#[:10]
            mask = ck_cd_exp_gene.isin(df['mapid'])
            existing_values = ck_cd_exp_gene[mask]
            print(existing_values)

        elif speci=='sorghum':
            input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/Sbi_TF_list.txt'  # 替换为真实的 csv 文件路径
            df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
            mapdata = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/sorghum/sorghum_geneid_trans.csv')
            gene_id_map = pd.Series(mapdata['merged_value'].values, index=mapdata['new_column'].values).to_dict()
            df['Gene_ID'] = df['Gene_ID'].map(gene_id_map)
            #ic(df)
            # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
            df['mapid'] = df.iloc[:,1]
            df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_mulgeneid/{data_type}_sorghum_tf_list.csv',index=False)
            
            ck_ce_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_sorghum/CK_CE/all_sorghum_CK_CE_DEG.csv')
            ck_ce_degree_gene = ck_ce_degree['Unnamed: 0']#[:10]
            mask = ck_ce_degree_gene.isin(df['mapid'])
            existing_values = ck_ce_degree_gene[mask]
            print(existing_values)

            ck_cd_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_sorghum/CK_CD/all_sorghum_CK_CD_DEG.csv')
            ck_cd_degree_gene = ck_cd_degree['Unnamed: 0']#[:10]
            mask = ck_cd_degree_gene.isin(df['mapid'])
            existing_values = ck_cd_degree_gene[mask]
            print(existing_values)

            ck_ce_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_sorghum/CK_CE/sorghum_CK_CE_DEG.csv')
            ck_ce_exp_gene = ck_ce_exp['FEATURE_ID']#[:10]
            mask = ck_ce_exp_gene.isin(df['mapid'])
            existing_values = ck_ce_exp_gene[mask]
            print(existing_values)

            ck_cd_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_sorghum/CK_CD/sorghum_CK_CD_DEG.csv')
            ck_cd_exp_gene = ck_cd_exp['FEATURE_ID']#[:10]
            mask = ck_cd_exp_gene.isin(df['mapid'])
            existing_values = ck_cd_exp_gene[mask]
            print(existing_values)

        else:
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
            #ic(relation)
            # 读取目标 CSV 文件中的某一列为 DataFrame
            input_csv_path = '/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/initial_data/indica/LOC_OS_TF_list.csv'  # 替换为真实的 csv 文件路径
            df =pd.read_csv(input_csv_path,sep='\t')  # 替换 ‘某一列’ 为你要读取的列名
            #ic(df)
            # 创建一个新的 DataFrame，包含 ID 和对应的 Relation
            df['mapid'] = df.iloc[:,1].map(relation).fillna("None")
            df.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/tf_mulgeneid/{data_type}_{speci}_tf_list.csv',index=False)

            ck_ce_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CE/all_{speci}_CK_CE_DEG.csv')
            ck_ce_degree_gene = ck_ce_degree['Unnamed: 0'][:10]
            mask = ck_ce_degree_gene.isin(df['mapid'])
            existing_values = ck_ce_degree_gene[mask]
            print(existing_values)

            ck_cd_degree=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/GRN_degree_DEG/{data_type}_{speci}/CK_CD/all_{speci}_CK_CD_DEG.csv')
            ck_cd_degree_gene = ck_cd_degree['Unnamed: 0'][:10]
            mask = ck_cd_degree_gene.isin(df['mapid'])
            existing_values = ck_cd_degree_gene[mask]
            print(existing_values)

            ck_ce_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/{speci}_CK_CE_DEG.csv')
            ck_ce_exp_gene = ck_ce_exp['FEATURE_ID'][:10]
            mask = ck_ce_exp_gene.isin(df['mapid'])
            existing_values = ck_ce_exp_gene[mask]
            print(existing_values)

            ck_cd_exp=pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/{speci}_CK_CD_DEG.csv')
            ck_cd_exp_gene = ck_cd_exp['FEATURE_ID'][:10]
            mask = ck_cd_exp_gene.isin(df['mapid'])
            existing_values = ck_cd_exp_gene[mask]
            print(existing_values)