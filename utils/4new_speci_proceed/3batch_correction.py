import pandas as pd
import os
from combat.pycombat import pycombat
from icecream import install,ic
import matplotlib.pyplot as plt
install()
import numpy as np

#
species=['wheat']#sorghum'','indica_BGI','indica','japonica',
for speci in species:
    folder_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/expression/{speci}'  # 设置你的文件夹路径

    # 读取所有文件
    all_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.txt')]
    #batch_labels = []  # 存储每个样本的批次信息
    data_frames = []  # 存储每个文件的DataFrame

    for idx, file_path in enumerate(all_files):
        df = pd.read_csv(file_path, sep='\t')  # 假设每个文件有 header
        #ic(df)
        df.set_index('FEATURE_ID', inplace=True)
        df=df+1e-16
        #df=np.log1p(df)
        data_frames.append(df)
        # batch_labels += [idx] * len(df)  # 假设每个文件的所有样本属于同一个批次


    # 合并所有数据框
    combined_df = pd.concat(data_frames, axis=1)
    ic(combined_df)
    threshold = 2/3

    row_threshold = combined_df.shape[1] * threshold
    combined_df = combined_df[combined_df.apply(lambda row: (row <= 1).sum() <= row_threshold, axis=1)]
    ic(combined_df)
    

    # # 删除超过2/3的数小于等于1的列
    # col_threshold = combined_df.shape[0] * threshold
    # combined_df = combined_df.loc[:, combined_df.apply(lambda col: (col <= 1).sum() <= col_threshold, axis=0)]
    
    # 筛选不满足删除条件的行（即小于等于1的元素个数不超过阈值的行）
    threshold_ratio = 2 / 3
    row_threshold = combined_df.shape[1] * threshold_ratio
    combined_df = combined_df.loc[(combined_df <= 1).sum(axis=1) <= row_threshold]

    # plt.boxplot(combined_df)
    # plt.savefig('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/utils/batchcorrection_former.png',bbox_inches='tight')
    ic(combined_df)
    combined_df = combined_df.dropna()
    # 进行批次效应校正

    batch = []
    for j in range(len(data_frames)):
        batch.extend([j for _ in range(len(data_frames[j].columns))])

    # run pyComBat
    df_corrected = pycombat(combined_df,batch)

    #df_corrected[df_corrected < 0] = 1e-10
    ic(df_corrected)
    # plt.boxplot(df_corrected)
    # plt.savefig('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/utils/batchcorrection_after.png',bbox_inches='tight')
    #ic(corrected_df)
    df_corrected.to_pickle(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression/{speci}_expression_batchcorrection.pkl')

#高粱数据由于各个批次差异太大，导致标准化后的数据无差异

#===============================================================================================================

species=['wheat']#'sorghum','indica_BGI','indica','japonica',
for speci in species:
    folder_path = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/count/{speci}'  # 设置你的文件夹路径

    # 读取所有文件
    all_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.txt')]
    #batch_labels = []  # 存储每个样本的批次信息
    data_frames = []  # 存储每个文件的DataFrame

    for idx, file_path in enumerate(all_files):
        df = pd.read_csv(file_path, sep='\t')  # 假设每个文件有 header
        #ic(df)
        df.set_index('FEATURE_ID', inplace=True)
        df=df+1e-16
        #df=np.log1p(df)
        data_frames.append(df)
        # batch_labels += [idx] * len(df)  # 假设每个文件的所有样本属于同一个批次


    # 合并所有数据框
    combined_df = pd.concat(data_frames, axis=1)
    ic(combined_df)
    threshold = 2/3

    row_threshold = combined_df.shape[1] * threshold
    combined_df = combined_df[combined_df.apply(lambda row: (row <= 1).sum() <= row_threshold, axis=1)]
    ic(combined_df)
    

    # # 删除超过2/3的数小于等于1的列
    # col_threshold = combined_df.shape[0] * threshold
    # combined_df = combined_df.loc[:, combined_df.apply(lambda col: (col <= 1).sum() <= col_threshold, axis=0)]
    
    # 筛选不满足删除条件的行（即小于等于1的元素个数不超过阈值的行）
    threshold_ratio = 2 / 3
    row_threshold = combined_df.shape[1] * threshold_ratio
    combined_df = combined_df.loc[(combined_df <= 1).sum(axis=1) <= row_threshold]

    # plt.boxplot(combined_df)
    # plt.savefig('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/utils/batchcorrection_former.png',bbox_inches='tight')
    ic(combined_df)
    combined_df = combined_df.dropna()
    # 进行批次效应校正

    batch = []
    for j in range(len(data_frames)):
        batch.extend([j for _ in range(len(data_frames[j].columns))])

    # run pyComBat
    df_corrected = pycombat(combined_df,batch)

    #df_corrected[df_corrected < 0] = 1e-10
    ic(df_corrected)
    # plt.boxplot(df_corrected)
    # plt.savefig('/home/win/4T/GeneDL/OSDrought_GCNAT_Link/utils/batchcorrection_after.png',bbox_inches='tight')
    #ic(corrected_df)
    df_corrected.to_pickle(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/count/{speci}_count_batchcorrection.pkl')
