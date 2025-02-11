import omicverse as ov
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from icecream import install,ic
install()
import pickle
import os
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white')

data_types=['TF','PPI','KEGG']
for data_type in data_types:
    species=['homo']
    #for organ in organs:
    for speci in species:
        test=0
        basepath='/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/'
        #gene_file=pd.read_csv(("/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/multispecies/TF_regulation_filter_data/zeamays/ALL_Genefile.csv"))
        #ov.utils.download_geneid_annotation_pair()
        #data=pd.read_csv('https://raw.githubusercontent.com/Starlitnightly/omicverse/master/sample/counts.txt',index_col=0,sep='\t',header=1)
        ref_data=pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/japonica_1134/RNASeq/DEG_data.csv')  #count数据
        ic(ref_data)
        if test==1:
            data=ref_data
        else:
            with open(f'{basepath}/{speci}_label_sample_count.pkl','rb') as fp:
                data=pickle.load(fp)
            if speci != 'homo':
                data.reset_index(inplace=True)
            new_columns = data.iloc[0]
            data.columns = new_columns
            data = data.iloc[1:].copy()  # 跳过第一行数据
            #ic(data)
            data = data.fillna(0)
            #if speci=='homo':
            new_columns = data.columns.tolist()
            new_columns[0] = 'FEATURE_ID'
            data.columns = new_columns
            data.iloc[:, 1:] = data.iloc[:, 1:].to_numpy().astype(int)
            #else:
                
                #data = data.astype({col: int for col in data.columns[1:]})  # 除了第一列外的列转换为int类型
            # 重置索引，这里直接在构造新DataFrame时避免了不必要的操作，因为新DataFrame默认已有连续的索引
            data.reset_index(drop=True, inplace=True)

            counter = {}
            new_columns = []
            for col in data.columns:
                counter[col] = counter.get(col, 0) + 1
                new_col_name = f"{col}.{counter[col]}" if counter[col] > 1 else col
                new_columns.append(new_col_name)
            #data.columns = new_columns
            #ic(new_columns)
            data=pd.DataFrame(data.values,index=list(data.index),columns=new_columns)
        ic(data.shape)

        data.iloc[:,1:] = data.iloc[:,1:].astype('float') + 1e-16

        ic(data)

        #data = data[data['FEATURE_ID'].isin(gene_file['Gene'])]
        ic(data.shape)

        '''示例数据格式  count数据
        GeneID,CD,CD,CD,CD,CD
        Os01g0100100,521,762,833,672
        Os01g0100200,1,0,0,2,0,0,0
        '''

        cd_columns = data.filter(like='CD', axis=1)
        ck_columns = data.filter(like='CK', axis=1)
        ce_columns = data.filter(like='CE', axis=1)
        ic(ck_columns)
        if test==0:
            cd_columns.reset_index(inplace=True)
            ck_columns.reset_index(inplace=True)
            #ce_columns.reset_index(inplace=True)
            cd_columns.drop(columns=cd_columns.columns[0], inplace=True)
            ck_columns.drop(columns=ck_columns.columns[0], inplace=True)
            #ce_columns.drop(columns=ce_columns.columns[0], inplace=True)

        #ic(ck_columns,cd_columns)
        #ic(cd_columns.shape,ck_columns.shape,ce_columns.shape)
        deg_data_ck_cd =  pd.concat([ck_columns, cd_columns], axis=1).astype(float)
        deg_data_ck_ce =  pd.concat([ck_columns, ce_columns], axis=1).astype(float)
        deg_data_cd_ce = pd.concat([cd_columns, ce_columns], axis=1).astype(float)

        dds_ck_cd=ov.bulk.pyDEG(deg_data_ck_cd)
        dds_ck_cd.drop_duplicates_index()
        print('... drop_duplicates_index success')
        dds_ck_cd.normalize()
        print('... estimateSizeFactors and normalize success')

        dds_ck_ce=ov.bulk.pyDEG(deg_data_ck_ce)
        dds_ck_ce.drop_duplicates_index()
        print('... drop_duplicates_index success')
        dds_ck_ce.normalize()
        print('... estimateSizeFactors and normalize success')

        dds_cd_ce=ov.bulk.pyDEG(deg_data_cd_ce)
        dds_cd_ce.drop_duplicates_index()
        print('... drop_duplicates_index success')
        dds_cd_ce.normalize()
        print('... estimateSizeFactors and normalize success')

        # cd_columns = cd_columns.astype('float16')
        # ce_columns = ce_columns.astype('float16')
        # ck_columns = ck_columns.astype('float16')
        # ic(cd_columns,ck_columns,ce_columns)
        treatment_groups_cd=cd_columns.columns#.tolist()
        treatment_groups_ce=ce_columns.columns#.tolist()
        control_groups_ck=ck_columns.columns#.tolist()
        #ic(treatment_groups_cd,treatment_groups_ce,control_groups_ck)
        #ic(dds_ck_cd)
        result_ck_cd=dds_ck_cd.deg_analysis(treatment_groups_cd,control_groups_ck,method='ttest')
        #ic(result_ck_cd)
        result_ck_cd = pd.DataFrame(result_ck_cd)
        combined_df_ck_cd = pd.concat([data.iloc[:, [0]], result_ck_cd], axis=1)

        result_ck_ce=dds_ck_ce.deg_analysis(treatment_groups_ce,control_groups_ck,method='ttest')
        result_ck_ce = pd.DataFrame(result_ck_ce)
        combined_df_ck_ce = pd.concat([ data.iloc[:, [0]], result_ck_ce], axis=1)

        result_cd_ce=dds_cd_ce.deg_analysis(treatment_groups_cd,treatment_groups_ce,method='ttest')
        result_cd_ce = pd.DataFrame(result_cd_ce)
        combined_df_cd_ce = pd.concat([ data.iloc[:, [0]], result_cd_ce], axis=1)

        def label_data(row):
            if row['sig'] == 'sig' and row['log2FC'] > 0:
                return 1
            elif row['sig'] == 'sig' and row['log2FC'] < 0:
                return -1
            else:
                return 0

        directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CD/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        combined_df_ck_cd.to_csv(f'{directory}/{speci}_CK_CD_DEG_alloutput.csv',index=False)
        ic(combined_df_ck_cd)
        df = combined_df_ck_cd
        #df['-log10pvalue'] = -np.log10(df['pvalue'])
        # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
        #selected_genes = df[(df['-log10pvalue'] > 0.05) & (df['log2FC'].abs() > 0.585)]
        # selected_genes_up = df[(df['pvalue'] < 0.05) & (df['abs(log2FC)'] > 0.585)]  #1.5倍
        # sorted_selected_genes_up = selected_genes_up.sort_values(by='abs(log2FC)', ascending=False)
        # sorted_selected_genes_up.to_csv(f'{directory}/{go}_{speci}_CK_CD_DEG.csv',index=False)

        selected_genes_up = df[(df['pvalue'] < 0.05) & (df['log2FC'] > 0.379)]
        selected_genes_up.to_csv(f'{directory}/{speci}_CK_CD_DEG_up.csv',index=False)
        selected_genes_down = df[(df['pvalue'] < 0.05) & (df['log2FC'] < -0.379)]
        selected_genes_down.to_csv(f'{directory}/{speci}_CK_CD_DEG_down.csv',index=False)

        all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
        sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
        sorted_selected_genes_up.to_csv(f'{directory}/{speci}_CK_CD_DEG.csv',index=False)

        #打印结果
        print(f"{speci}CK_CD上调基因数量：{len(selected_genes_up)}")
        print(f"{speci}CK_CD下调基因数量：{len(selected_genes_down)}")

        directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CK_CE/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        combined_df_ck_ce['Label'] = combined_df_ck_ce.apply(label_data, axis=1)
        # exp_data = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv')
        # exp_data['Label'] =combined_df_ck_ce.apply(label_data, axis=1)
        combined_df_ck_ce.to_csv(f'{directory}/{speci}_CK_CE_DEG_alloutput.csv',index=False)

        df = combined_df_ck_ce
        df['-log10pvalue'] = -np.log10(df['pvalue'])
        # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
        selected_genes_up = df[(df['pvalue'] < 0.05) & (df['log2FC'] > 0.379)]
        selected_genes_up.to_csv(f'{directory}/{speci}_CK_CE_DEG_up.csv',index=False)
        selected_genes_down = df[(df['pvalue'] < 0.05) & (df['log2FC'] < -0.379)]
        selected_genes_down.to_csv(f'{directory}/{speci}_CK_CE_DEG_down.csv',index=False)
        print(f"{speci}CK_CE上调基因数量：{len(selected_genes_up)}")
        print(f"{speci}CK_CE下调基因数量：{len(selected_genes_down)}")
        all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
        sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
        sorted_selected_genes_up.to_csv(f'{directory}/{speci}_CK_CE_DEG.csv',index=False)

        directory=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/model_result/EXP_DEG/{data_type}_{speci}/CD_CE/'
        if not os.path.exists(directory):
            os.makedirs(directory)
        combined_df_cd_ce['Label'] = combined_df_cd_ce.apply(label_data, axis=1)
        # exp_data = pd.read_csv('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/data/RNASeq/Train_validation_test/ExpressionData_unique_networkfilter.csv')
        # exp_data['Label'] =combined_df_ck_ce.apply(label_data, axis=1)
        combined_df_cd_ce.to_csv(f'{directory}/{speci}_CD_CE_DEG_alloutput.csv',index=False)

        df = combined_df_cd_ce
        df['-log10pvalue'] = -np.log10(df['pvalue'])
        # 筛选-log10pvalue大于0.05的行，并进一步筛选出abs(log2FC) > 0.585的行
        #selected_genes = df[(df['-log10pvalue'] > 0.05) & (df['log2FC'].abs() > 0.585)]
        # 根据abs(log2FC)降序排列
        selected_genes_up = df[(df['pvalue'] < 0.05) & (df['log2FC'] > 0.379)]
        selected_genes_up.to_csv(f'{directory}/{speci}_CD_CE_DEG_up.csv',index=False)
        selected_genes_down = df[(df['pvalue'] < 0.05) & (df['log2FC'] < -0.379)]
        selected_genes_down.to_csv(f'{directory}/{speci}_CD_CE_DEG_down.csv',index=False)
        print(f"{speci}CD_CE上调基因数量：{len(selected_genes_up)}")
        print(f"{speci}CD_CE下调基因数量：{len(selected_genes_down)}")
        all_select_gene_up = pd.concat([selected_genes_up,selected_genes_down],axis=0)
        sorted_selected_genes_up = all_select_gene_up.sort_values(by='abs(log2FC)', ascending=False)
        sorted_selected_genes_up.to_csv(f'{directory}/{speci}_CD_CE_DEG.csv',index=False)




