import numpy as np
from icecream import install,ic
install()



# 替换下面的路径为你的 .npy 文件的实际路径
H_T_bulk = histogram2d(x_tf_bulk, x_gene_bulk, bins=32)

file_path = '/mnt/h/study/deep_learning/OSDrought_GCNAT_Link/data/NEPDF/img/zdata_tf2.npy'

# 加载数据
data = np.load(file_path)

# 打印数据
ic(data)