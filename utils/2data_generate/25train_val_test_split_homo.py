#train 96459  val 14401  test  33230
#7：1：2

import pandas as pd
from sklearn.model_selection import train_test_split

# # 读取CSV文件
# df = pd.read_csv('/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/all_labels_unique.csv')

# # 划分正负样本
# positive_samples = df[df['Label'] == 1]
# negative_samples = df[df['Label'] == 0]

# # 划分正样本为训练集、验证集和测试集
# X_pos_train, X_pos_remain, y_pos_train, y_pos_remain = train_test_split(
#     positive_samples, positive_samples['Label'], test_size=0.3, stratify=positive_samples['Label'], random_state=42)
# X_pos_val, X_pos_test, y_pos_val, y_pos_test = train_test_split(
#     X_pos_remain, y_pos_remain, test_size=0.2, stratify=y_pos_remain, random_state=42)

# # 划分负样本为训练集、验证集和测试集
# X_neg_train, X_neg_val, y_neg_train, y_neg_val = train_test_split(
#     negative_samples, negative_samples['Label'], test_size=0.3, stratify=negative_samples['Label'], random_state=42)
# X_neg_val, X_neg_test, y_neg_val, y_neg_test = train_test_split(
#     X_neg_val, y_neg_val, test_size=0.2, stratify=y_neg_val, random_state=42)

# # 合并正负样本的训练集、验证集和测试集
# X_train = pd.concat([X_pos_train, X_neg_train], ignore_index=True)
# y_train = pd.concat([y_pos_train, y_neg_train], ignore_index=True)
# X_val = pd.concat([X_pos_val, X_neg_val], ignore_index=True)
# y_val = pd.concat([y_pos_val, y_neg_val], ignore_index=True)
# X_test = pd.concat([X_pos_test, X_neg_test], ignore_index=True)
# y_test = pd.concat([y_pos_test, y_neg_test], ignore_index=True)

# # 打印各数据集中正负样本的数量
# print("训练集正样本数量:", len(X_train[y_train == 1]))
# print("训练集负样本数量:", len(X_train[y_train == 0]))
# print("验证集正样本数量:", len(X_val[y_val == 1]))
# print("验证集负样本数量:", len(X_val[y_val == 0]))
# print("测试集正样本数量:", len(X_test[y_test == 1]))
# print("测试集负样本数量:", len(X_test[y_test == 0]))

# # 保存划分后的数据集为CSV文件
# X_train.to_csv('/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/Train_validation_test/Train_set.csv')
# X_val.to_csv('/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/Train_validation_test/Validation_set.csv')
# X_test.to_csv('/mnt/h/study/deep_learning/gene/project/GNNLink/OsData/Train_validation_test/Test_set.csv')

#----------------------------------------------第二种划分方式-------------------------------------------------
import pandas as pd
from sklearn.model_selection import train_test_split
from icecream import install,ic
install()
# 读取原始数据集的 CSV 文件'indica',
for speci in ['TF','PPI','KEGG']:
    datapath=f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{speci}/homo/'
    df = pd.read_csv(f'{datapath}/all_labels.csv')

    ic(len(df))
    # 等间隔取样得到测试集
    trainval_df, test_df = train_test_split(df, test_size=0.2, random_state=42)

    # test_interval = int(len(df)/7000)
    # test_indices = list(range(0, len(df), test_interval))
    # test_df = df.iloc[test_indices]

    # # 从原始数据集中删除测试集的行，剩下的是训练验证集
    # trainval_df = df.drop(test_indices)

    # 随机划分训练集和验证集
    train_df, val_df = train_test_split(trainval_df, test_size=0.25, random_state=42)

    # 打印划分后的数据集大小
    print('训练集大小:', len(train_df))
    print('验证集大小:', len(val_df))
    print('测试集大小:', len(test_df))

    train_neg_samples = train_df[train_df['Label'] == 0].shape[0]
    train_pos_samples = train_df[train_df['Label'] == 1].shape[0]

    # 统计验证集中标签为0和1的样本个数
    val_neg_samples = val_df[val_df['Label'] == 0].shape[0]
    val_pos_samples = val_df[val_df['Label'] == 1].shape[0]

    # 统计测试集中标签为0和1的样本个数
    test_neg_samples = test_df[test_df['Label'] == 0].shape[0]
    test_pos_samples = test_df[test_df['Label'] == 1].shape[0]

    # 打印结果
    print("训练集负样本个数:", train_neg_samples)
    print("训练集正样本个数:", train_pos_samples)
    print("验证集负样本个数:", val_neg_samples)
    print("验证集正样本个数:", val_pos_samples)
    print("测试集负样本个数:", test_neg_samples)
    print("测试集正样本个数:", test_pos_samples)


    # 将划分后的数据集保存到不同的 CSV 文件
    train_df.to_csv(f'{datapath}/Train_set.csv')
    val_df.to_csv(f'{datapath}/Validation_set.csv')
    test_df.to_csv(f'{datapath}/Test_set.csv')























#测试集手动划分，隔n个取一个，再划分训练和测试
#检查正样本和负样本之间的数据有没有重合    有重合！！！
#np.ones    已改回



