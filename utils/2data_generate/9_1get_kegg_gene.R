library(curl)

library("XML")

library("methods")

library(KEGGgraph)
# KEGG的物种界面 https://www.genome.jp/kegg/tables/br08606.html
#https://www.jianshu.com/p/e6e2126b2001
#osa	KGB	Oryza sativa japonica (Japanese rice)
#zma	KGB	Zea mays (maize)

library(KEGGREST)
#=========================================下载kgml文件==================================================================
# # 定义物种代码，例如 "hsa" 代表人类
species_code <- "osa"

# 获取物种的所有路径
pathways <- keggList("pathway", species_code)
print(pathways)

# pathways_df <- data.frame(
#   ID = names(pathways),
#   Description = pathways,
#   stringsAsFactors = FALSE
# )
# file_path <- "/home/win/4T/GeneDL/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/japonica/pathways_list.csv"  # 替换为实际的文件路径
# write.csv(pathways_df, file = file_path, row.names = FALSE)

# 提取路径 ID
pathway_ids <- sub(paste0(species_code, ":"), "", names(pathways))
print(pathway_ids)
# 创建保存 KGML 文件的目录
output_dir <- paste0('/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/', "zma_kgml")
dir.create(output_dir, showWarnings = FALSE)

download_paths <- vector("character", length(pathway_ids))
# 下载所有路径的 KGML 文件
for (i in seq_along(pathway_ids)) {
  pathway_id <- pathway_ids[i]
  kgml_url <- paste0("http://rest.kegg.jp/get/", pathway_id, "/kgml")
  kgml_file <- paste0(output_dir, "/", pathway_id, ".xml")
  
  # 下载并保存 KGML 文件
  tryCatch({
    download.file(kgml_url, kgml_file, mode = "wb")
    # 保存下载路径
    download_paths[i] <- kgml_file
    # 打印下载状态
    cat("Downloaded:", kgml_file, "\n")
  }, error = function(e) {
    warning(sprintf("Failed to download %s: %s", kgml_url, e$message))
  })
}

cat("All KGML files have been downloaded to:", output_dir, "\n")
#=======================================得到基因名称================================================

# process_kgml_file <- function(input_file, output_file) {
#   # 解析 KGML 文件
#   mapkG <- parseKGML2Graph(input_file, expandGenes = TRUE, genesOnly = TRUE)
  
#   # 提取节点和边
#   mapkNodes <- nodes(mapkG)
#   mapkEdges <- edges(mapkG)
  
#   # 过滤空边
#   mapkEdges <- mapkEdges[sapply(mapkEdges, length) > 0]
  
#   # 如果 mapkEdges 为空，则跳过处理
#   if (length(mapkEdges) == 0) {
#     cat("No edges found in:", input_file, "\n")
#     return(NULL)
#   }
  
#   # 处理边信息
#   res <- lapply(1:length(mapkEdges), function(t) {
#     name <- names(mapkEdges)[t]
#     len <- length(mapkEdges[[t]])
#     do.call(rbind, lapply(1:len, function(n) {
#       c(name, mapkEdges[[t]][n])
#     }))
#   })
  
#   result <- data.frame(do.call(rbind, res))
  
#   # 保存结果到文件
#   write.table(result, output_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
# }

# # 定义遍历文件夹并处理文件的函数
# process_kgml_files_in_directory <- function(input_dir, output_dir) {
#   # 获取所有 xml 文件
#   kgml_files <- list.files(input_dir, pattern = "\\.xml$", full.names = TRUE)
  
#   # 创建输出目录（如果不存在）
#   if (!dir.exists(output_dir)) {
#     dir.create(output_dir, recursive = TRUE)
#   }
  
#   # 遍历所有文件并处理
#   for (file in kgml_files) {
#     # 构建输出文件路径
#     output_file <- file.path(output_dir, paste0(basename(file), "_edges.txt"))
    
#     # 处理文件并保存结果
#     process_kgml_file(file, output_file)
    
#     # 打印处理状态
#     cat("Processed and saved:", output_file, "\n")
#   }
# }

# # 指定输入和输出目录
# input_directory <- "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zma_kgml"
# output_directory <- "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/zma_edge"

# # 处理所有文件
# process_kgml_files_in_directory(input_directory, output_directory)
#==================================================================================================
# write.table(mapkNodes, "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/nodes.txt", sep = "\t", row.names = F, col.names = F, quote = F ,append = T)


#interested_pathway <- c('hsa04390','hsa04371')

# for (i in interested_pathway){

# url <- paste0("http://rest.kegg.jp/get/",i,"/kgml")

# assign(paste0(i,"_kgml"), tempfile())

# curl_download(url, paste0(i,"_kgml"))

# }

##############生成进行网络图绘制的边文件和节点文件#####################

# for (i in interested_pathway){

#  mapkG <- parseKGML2Graph(paste0(i,"_kgml"),expandGenes=TRUE, genesOnly = TRUE)

#  mapkNodes <- nodes(mapkG)

#  mapkEdges <- edges(mapkG)

#  mapkEdges <- mapkEdges[sapply(mapkEdges, length) > 0]

#  res <- lapply(1:length(mapkEdges), function(t){

#   name <- names(mapkEdges)[t]

#   len <- length(mapkEdges[[t]])

#   do.call(rbind, lapply(1:len, function(n){

#    c(name, mapkEdges[[t]][n])

#   }))

#  })

#  result <- data.frame(do.call(rbind, res))

#  write.table(result, "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/edges.txt", sep = "\t", row.names = F, col.names = F, quote = F ,append = T)

#  write.table(mapkNodes, "/home/win/16t1/Deeplearning/data_output/GRN_result/Batch_correction/data/batch_correction_data/KEGG/nodes.txt", sep = "\t", row.names = F, col.names = F, quote = F ,append = T)

# }










