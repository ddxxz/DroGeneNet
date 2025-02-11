import pandas as pd
import numpy as np
import h5py
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import sys
import pickle
# 生成二维数组的函数，这里只是创建随机数组作为演示
def create_data(idx,indexed_TFs,train_expression,indexed_Targets):
    gene_id = indexed_TFs.iloc[idx]  # idx 是你想要检索的行号
    #print(train_expression.shape)
    train_expression = train_expression.drop_duplicates(subset='FEATURE_ID')
    TF_train_data= train_expression.loc[train_expression['FEATURE_ID'] == gene_id]
    gene_id = indexed_Targets.iloc[idx]
    Target_train_data = train_expression.loc[train_expression['FEATURE_ID'] == gene_id]
    #ic(selected_train_data.shape)
    #print(np.array(Target_train_data).shape)
    data = np.histogram2d(np.array(TF_train_data).squeeze(axis=0)[1:], np.array(Target_train_data).squeeze(axis=0)[1:], bins=32)
    data = data[0].T

    return (idx, data)

def main():

    dataset = sys.argv[1]
    species = sys.argv[2]
    data_type = sys.argv[3]
    dataset_dir = f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{species}'

    train_data = pd.read_csv(f'{dataset_dir}/{dataset}_set.csv',index_col=0).values
    with open(f'{dataset_dir}/ExpressionData_unique_networkfilter.pkl', 'rb') as file:
    # 使用 pickle 加载数据
        train_expression = pickle.load(file)
    #train_expression = pd.read_csv(f'/home/win10/OSDrought_GCNAT_Link/data/{folder}/Train_validation_test/ExpressionData_unique_networkfilter.csv')
    gene_file = pd.read_csv(f'{dataset_dir}/ALL_Genefile.csv')

    train_input = train_data[:,:2]
    indexed_TFs = gene_file.set_index('index').reindex(train_input[:,0])['Gene']
    indexed_Targets = gene_file.set_index('index').reindex(train_input[:,1])['Gene']
    print("start.......")
    # 使用多进程池来加速数组生成
    with Pool(processes=10) as pool:
        # 使用functools.partial来给函数固定参数
        # 这里以10x10的数组为例
        create_data_partial = partial(create_data,indexed_TFs=indexed_TFs,train_expression=train_expression,indexed_Targets=indexed_Targets)
        
        # map将索引映射到函数上，并并行执行
        #all_data =     
    # # 保证生成数据的顺序与索引相符
    # all_data.sort(key=lambda x: x[0])

    # 将生成的数据写入HDF5文件
        with h5py.File(f'{dataset_dir}/NEPDF/{dataset}_data.h5', 'w') as hdf_file:
            print(len(train_data))
            for idx, data in tqdm(pool.map(create_data_partial, range(len(train_data)))):
                dataset_name = f'dset_{idx}'
                hdf_file.create_dataset(dataset_name, data=data)

if __name__ == '__main__':
    main()




