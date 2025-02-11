#install.packages("WGCNA")
#BiocManager::install('WGCNA')

#参考代码 https://blog.csdn.net/kanghua_du/article/details/129232316
#https://www.jianshu.com/p/cf7bd18612a3
#https://mp.weixin.qq.com/s/M0LAlE-61f2ZfpMiWN-iQg
#导出网络https://www.jianshu.com/p/cb83950d5deb

library(WGCNA)
options(stringsAsFactors = FALSE)
## 打开多线程
enableWGCNAThreads()
speci_list <- c('indica_BGI', 'japonica', 'zeamays')
treat_list <- c('CK', 'CD')

for (speci in speci_list) {
  for (treat in treat_list) {
    # 构建文件路径

    file_path <- sprintf("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_out/expression_WGCNA/%s_%s_WGCNA.txt", speci, treat)
    
    exr1_symbol_no_dup <- read.csv(file_path,sep='\t',row.names = 1)
    dim(exr1_symbol_no_dup)
    head(exr1_symbol_no_dup)
    colnames(exr1_symbol_no_dup)

    Rfile_path <- sprintf("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/expression_WGCNA/%s_%s_WGCNA.RData", speci, treat)
    if (file.exists(Rfile_path)) {
      # 如果文件存在，跳过此次循环
      next
    }
    #转置
    mydata <- exr1_symbol_no_dup
    datExpr2 = data.frame(t(exr1_symbol_no_dup))
    colnames(datExpr2) <- rownames(mydata)
    rownames(datExpr2) <- colnames(mydata)
    head(datExpr2)
    dim(datExpr2)
    datExpr <- datExpr2

    #========================================使用=================================================
    # library(extrafont)
    # font_import()  # 这一步可能需要一些时间，请耐心等待
    # loadfonts(device = "pdf")  # 或者 "win" 如果你在 Windows 上


    softPower = 6
    adjacency = adjacency(datExpr, power = softPower)
    # Turn adjacency into topological overlap
    TOM = TOMsimilarity(adjacency)
    dissTOM = 1-TOM 
    # Call the hierarchical clustering function 
    geneTree = hclust(as.dist(dissTOM), method = "average")
    # Plot the resulting clustering tree (dendrogram) 
    #sizeGrWindow(12,9) 
    #plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)

    # We like large modules, so we set the minimum module size relatively high: 
    minModuleSize = 30
    # Module identification using **dynamic tree cut**: 
    dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
    table(dynamicMods)
    # Convert numeric lables into colors 
    dynamicColors = labels2colors(dynamicMods)
    table(dynamicColors)
    # Plot the dendrogram and colors underneath
    #sizeGrWindow(8,6)
    #plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
    # Calculate eigengenes 
    MEList = moduleEigengenes(datExpr, colors = dynamicColors)
    MEs = MEList$eigengenes 
    # Calculate dissimilarity of module eigengenes
    MEDiss = 1-cor(MEs)
    # Cluster module eigengenes 
    METree = hclust(as.dist(MEDiss), method = "average")
    # Plot the result 
    sizeGrWindow(7, 6) 
    plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
    MEDissThres = 0.25 
    # Plot the cut line into the dendrogram 
    abline(h=MEDissThres, col = "red")
    # Call an automatic merging function 
    merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
    # The merged module colors 
    mergedColors = merge$colors 
    # Eigengenes of the new merged modules: 
    mergedMEs = merge$newMEs
    # sizeGrWindow(12, 9) 
    # #pdf(file = "Plots/geneDendro-3.pdf", wi = 9, he = 6)
    # png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/5_Dynamic Tree Cut.png", width = 500, height = 300)  # 640x480 pixels
    # par(cex.lab=2, cex.axis=2, cex.main=2, font.lab=2, font.axis=2)  #font.lab: 用于坐标轴标签的字体。   font.axis: 用于坐标轴刻度的字体。  cex.lab: 坐标轴标签的字体大小。 cex.axis: 坐标轴刻度文字的字体大小。 
    # plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) 
    # dev.off()
    # Rename to moduleColors 
    moduleColors = mergedColors
    # Construct numerical labels corresponding to the colors
    colorOrder = c("grey", standardColors(50))
    moduleLabels = match(moduleColors, colorOrder)-1
    MEs = mergedMEs
    # Save module colors and labels for use in subsequent parts
    save(MEs, moduleLabels, moduleColors,dynamicColors,dynamicMods,geneTree,dissTOM,TOM,file = sprintf("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/expression_WGCNA/%s_%s_WGCNA.RData",speci,treat))

    load(sprintf("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/expression_WGCNA/%s_%s_WGCNA.RData",speci,treat))
    #par(family="Times New Roman", cex=2)
    png(sprintf("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/expression_WGCNA/%s_%s_WGCNA.png",speci,treat), width = 500, height = 300)  # 640x480 pixels
    par(cex.lab=2, cex.axis=2, cex.main=2, font.lab=2, font.axis=2)  #font.lab: 用于坐标轴标签的字体。   font.axis: 用于坐标轴刻度的字体。  cex.lab: 坐标轴标签的字体大小。 cex.axis: 坐标轴刻度文字的字体大小。 
    plotDendroAndColors(geneTree, cbind(dynamicColors, moduleColors), 
                        groupLabels = c("Unmerged", "Merged"),
                        cex.colorLabels=1.2,
                        #cex.rowText=2,
                        dendroLabels = FALSE, hang = 0.03, 
                        addGuide = TRUE, guideHang = 0.05)
    # mtext(side = 2, line = 2, at = -0.5, "Dynamic Tree Cut", cex = 2)
    # mtext(side = 2, line = 2, at = -1.5, "Merged dynamic", cex = 2)
    dev.off()

    #输出所有modules
    modules = unique(moduleColors)
    # 获取所有基因名称
    probes = names(datExpr)
    # 初始化一个空列表，用于存储所有模块的数据
    allModulesData <- list()
    threshold = 0.02
    # 遍历每个模块颜色
    for (module in modules) {
    # 选取module中的基因
    inModule = (moduleColors == module)
    modProbes = probes[inModule]
    # 选择对应的拓扑重叠矩阵
    modTOM = TOM[inModule, inModule]
    dimnames(modTOM) = list(modProbes, modProbes)
    # 转换为数据框并添加模块信息
    modDataFrame = as.data.frame(as.table(modTOM))
    modDataFrame = modDataFrame[modDataFrame$Freq > threshold, ]
    modDataFrame$Module = module
    # 将当前模块的数据添加到列表中
    allModulesData[[module]] <- modDataFrame
    }
    # 将所有模块的数据合并成一个数据框
    combinedData <- do.call(rbind, allModulesData)
    # 输出为CSV文件
    write.csv(combinedData, file = sprintf("/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/out/expression_WGCNA/%s_%s_WGCNA.csv",speci,treat), row.names = FALSE, quote = FALSE)
  }
}
#基因过滤
# datExpr1<-datExpr2
# gsg = goodSamplesGenes(datExpr1, verbose = 3);
# gsg$allOK
# if (!gsg$allOK){
#   # Optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0) 
#     printFlush(paste("Removing genes:", paste(names(datExpr1)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0) 
#     printFlush(paste("Removing samples:", paste(rownames(datExpr1)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   datExpr1 = datExpr1[gsg$goodSamples, gsg$goodGenes]
# }

# #绘制样本聚类图
# sampleTree = hclust(dist(datExpr1), method = "average")
# png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/1_sample_clustering.png", width = 480, height = 320)  # 480x320 pixels
# par(cex = 0.7)
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# dev.off()

# # Plot the sample clustering with outliers removed
# png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/2_sample_clustering_delete_outliers.png", width = 640, height = 480)  # 640x480 pixels
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# abline(h = 100000, col = "red")  # Modify the 'h' parameter as needed for your data
# dev.off()

# clust = cutreeStatic(sampleTree, cutHeight = 100000, ##cutHeight依据自己的数据设置 
#                      minSize = 10)
# keepSamples = (clust==1)
# datExpr = datExpr1[keepSamples, ]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# dim(datExpr)
# head(datExpr)
####
#datExpr0 <- datExpr
#此处可以和表型相关联 干旱设置为1，正常设置为0
############### 载入性状数据## input trait data###############
# traitData = read.csv("TraitData.csv",row.names=1)
# head(traitData)
# allTraits = traitData
# dim(allTraits)
# names(allTraits)
# # 形成一个类似于表达数据的数据框架
# fpkmSamples = rownames(datExpr0)
# traitSamples =rownames(allTraits)
# traitRows = match(fpkmSamples, traitSamples)
# datTraits = allTraits[traitRows,]
# rownames(datTraits)
# collectGarbage()
#=========二次样本聚类===============================
# sampleTree2 = hclust(dist(datExpr), method = "average")
# traitColors = numbers2colors(datTraits, signed = FALSE)
# pdf(file="3_Sample_dendrogram_and_trait_heatmap.pdf",width=8,height=6)
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits),
#                     main = "Sample dendrogram and trait heatmap",cex.colorLabels = 1.5, cex.dendroLabels = 1, cex.rowText = 2)
# dev.off()

#筛选软阈值
#enableWGCNAThreads()
# # Choose a set of soft-thresholding powers
# #powers = c(1:30)
# powers = c(c(1:10), seq(from = 12, to=30, by=2))
# # Call the network topology analysis function
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# #绘图soft power plot
# png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/4_软阈值选择.png", width = 640, height = 480)  # 640x480 pixels
# #pdf(file="4_软阈值选择.pdf",width=12,height= 8)
# par(mfrow = c(1,2))
# cex1 = 0.85
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# # this line corresponds to using an R^2 cut-off of h
# abline(h=0.85,col="red")
# # Mean connectivity as a function of the soft-thresholding power
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# dev.off()

# #softPower =sft$powerEstimate
# sft$powerEstimate
# softPower = 14

# # #模块可视化
# net = blockwiseModules(datExpr, power = 6,#手动改power
#                        #signed, unsigned
#                        TOMType = "unsigned", minModuleSize = 30,#20, 25 
#                        reassignThreshold = 0, mergeCutHeight = 0.25, #mergecutheight 0.25
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,maxBlockSize = 2000,
#                        saveTOMFileBase = "MyTOM",
#                        verbose = 3)



# color<-unique(moduleColors)
# for (i  in 1:length(color)) {
#   y=t(assign(paste(color[i],"expr",sep = "."),datExpr[moduleColors==color[i]]))
#   write.csv(y,paste('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/6',color[i],"csv",sep = "."),quote = F)
# }

# mergedData <- data.frame()
# # 遍历每个颜色模块，合并数据
# for (color in colors) {
#   moduleData <- datExpr[, moduleColors == color, drop = FALSE]
#   colnames(moduleData) <- paste(color, colnames(moduleData), sep = "_")
#   mergedData <- cbind(mergedData, moduleData)
# }
# # 写入CSV文件
# outputFilePath <- "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/merged_modules.csv"
# write.csv(mergedData, outputFilePath, quote = FALSE)

# table(net$colors)
# #绘制模块聚类图
# mergedColors = labels2colors(net$colors)
# table(mergedColors)
# png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/5_Dynamic Tree Cut.png", width = 640, height = 480)  # 640x480 pixels
# #pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6)
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()

# #模块的合并
# # 合并
# merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# # The merged module colors
# mergedColors = merge$colors
# # Eigengenes of the new merged modules:
# mergedMEs = merge$newMEs
# table(mergedColors)

# #sizeGrWindow(12, 9)
# png("/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/7_merged dynamic.png", width = 640, height = 480)  # 640x480 pixels
# #pdf(file="7_merged dynamic.pdf", width = 9, height = 6)
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
#                     c("Dynamic Tree Cut", "Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# dev.off()

#输出所有的模块基因
# moduleLabels = net$colors
# moduleColors = labels2colors(net$colors)
# MEs = net$MEs
# geneTree = net$dendrograms[[1]]
# #输出所有modules

# color<-unique(moduleColors)
# for (i  in 1:length(color)) {
#   y=t(assign(paste(color[i],"expr",sep = "."),datExpr[moduleColors==color[i]]))
#   write.csv(y,paste('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/6',color[i],"csv",sep = "."),quote = F)
# }

# mergedData <- data.frame()
# # 遍历每个颜色模块，合并数据
# for (color in colors) {
#   moduleData <- datExpr[, moduleColors == color, drop = FALSE]
#   colnames(moduleData) <- paste(color, colnames(moduleData), sep = "_")
#   mergedData <- cbind(mergedData, moduleData)
# }
# # 写入CSV文件
# outputFilePath <- "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/out/merged_modules.csv"
# write.csv(mergedData, outputFilePath, quote = FALSE)

# save.image(file = "module_splitted.RData")
# load("module_splitted.RData")

#模块和表型数据的相关性热图
## 表型
# #samples <- read.csv("TraitData.csv",row.names = 1,header = T)
# samples <- traitData
# samples <- samples[, -(6:6)]
# print(samples)
# ### ----------------------------------------------------------------------------
# ##        (最重要的) 模块和性状的关系
# moduleLabelsAutomatic <-  net$colors
# moduleColorsAutomatic <-  labels2colors(moduleLabelsAutomatic)
# moduleColorsWW <-  moduleColorsAutomatic
# MEs0 <-  moduleEigengenes(datExpr, moduleColorsWW)$eigengenes
# ## 赋值，后续可能用得到
# moduleColors = moduleColorsWW
# ####
# MEsWW <-  orderMEs(MEs0)
# modTraitCor <-  cor(MEsWW, samples, use = "p")
# colnames(MEsWW)
# ###赋值
# modlues = MEsWW
# #write.csv(modlues,file = "./modules_expr.csv")
# modTraitP <-  corPvalueStudent(modTraitCor, nSamples)
# textMatrix <- paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
# dim(textMatrix) <-  dim(modTraitCor)








