#!/bin/bash

#=== 环境 r_env ===
# 定义植物种类和处理类型
plants="indica_BGI japonica zeamays"
treatments="CK CD "  #CE

# 外层循环遍历植物种类
for plant in $plants; do
    for treatment in $treatments; do
        # 创建输出目录，使用mkdir -p可以创建多级目录，即使父目录不存在
        out_dir="/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/WGCNA_R/out/$plant/$treatment/"
        mkdir -p "$out_dir"
        
        echo "Running R script for $plant with $treatment at layer $layer"
        Rscript /home/win/4T/GeneDL/DXZ_DL/expriments/GRN/utils/data_proceed/9WGCNA.R \
                -i /home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/WGCNA_R/input/$plant/input_$treatment.json \
                -o "$out_dir"
    done
done

#Rscript /home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/src/WGCNA.R -i /home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_WGCNA/indica/input/input.json -o /home/win/4T/GeneDL/OSDrought_GCNAT_Link/plot/multispecies_WGCNA/indica/out/




















