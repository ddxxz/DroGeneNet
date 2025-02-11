import h5py
import pandas as pd
import numpy as np
from icecream import install,ic
install()
import sys
import pickle

data_type = sys.argv[1]
speci = sys.argv[2]
#-------------------------------------生成表达数据---------------------------------------------------------
with open(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/ExpressionData_unique_networkfilter.pkl', 'rb') as file:
    expression_data = pickle.load(file)
#expression_data = pd.read_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/ExpressionData_unique_networkfilter.csv')
expression_data = expression_data.T
ic(expression_data)

#print(expression_data)
#print(expression_data.iloc[0,:])
x= pd.DataFrame([i+1 for i in range(len(expression_data))])
#print(x)
#print(expression_data.columns)

axis0 = np.array(expression_data.iloc[0,:])
axis1 = pd.DataFrame(np.array([i+1 for i in range(len(expression_data))]).astype('int'))     
block0_items = pd.DataFrame(np.array(expression_data.iloc[0, :]).astype('S'))
block0_values = pd.DataFrame(np.array(expression_data.iloc[1:, :]).astype('float'))
ic(axis1)
new_data = pd.DataFrame(np.array(expression_data.iloc[1:, :]).astype('float'),columns=axis0)
ic(new_data)
new_data.to_csv(f'/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/{data_type}/{speci}/NEPDF/ExpressionData.csv',index=False)

