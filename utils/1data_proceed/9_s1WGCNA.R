#Rscript script_name.R arg1 arg2
# install.packages("devtools", repos='https://cran.r-project.org')
# devtools::install_version("ggplot2", version="2.2.1")
# install.packages("optparse", repos='https://cran.r-project.org')
# install.packages("jsonlite", repos='https://cran.r-project.org')
# install.packages('amap', repos='https://cran.r-project.org')
# install.packages("dendextend", repos='https://cran.r-project.org')
# install.packages("WGCNA", repos='https://cran.r-project.org')

library(optparse)
option_list <- list(
  make_option(c("-i", "--infile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/tool_runing/8ed926747e04570e4a5345bfa133185b/input.json',
              action = "store", help = "Input a exp file path!"
  ),
  make_option(c("-o", "--outfile"), type = "character", default = '/pub1/data/mg_projects/projects/web_script/web_file_catche/runing/dae42f7681eca28355c6ba9105a3b5ff',
              action = "store", help = "Input a outfolder path!"
  )
)
logs=c()

tryCatch({
  Args <- commandArgs()
  opt = parse_args(OptionParser(option_list = option_list, usage = "Data press"))
  #logs=c(logs,paste0('geting data:',paste0(paste0(names(opt),'=',opt),collapse = ',')))
  logs=c(logs,paste0('run wgcna.R-',basename(opt$outfile)))

  library(jsonlite)
  data<-jsonlite::stream_in(file(opt$infile),pagesize = 1000)
  
#exp_path='/pub1/data/user_data/13456826965/GSE1561/GSE1561_GPL96_sample_exp.txt'
#geneMode='mad'
#geneModeValue='50'
#sampleCut='1'
#netType='unsigned'
#geneModeValue=as.numeric(geneModeValue)
#outFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/WGCNA1'
#minModuleSize=30
#deepSplit=4
#mergeCutHeight=0.5

#oldFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/test_data/WGCNA'
method=unlist(data$method)#'cli'#module
#cliDatas=list()
#list(A=list(name:xxx,sample:c(),value:c()),B=list())
oldFolder=unlist(data$oldFolder)

library('amap')
library(dendextend)
getTreeNode=function(hc_obj){
  dend=as.dendrogram(hc_obj)
  capture.output(str(dend))->stx2
  #write.table(stx,row.names = F,col.names = F ,paste0(outFolder,'/wgcna_dissTOM_tree.text'))
  gsub('^ \\s+','',stx2)->stx2
  nodeMap=rbind()
  leafs=rbind()
  all_node_obj=rbind()
  for(i in length(stx2):1){
    #  i=length(stx2)
    #i=547
    if(length(grep('--leaf',stx2[i]))==0){
      h=unlist(stringr::str_split(stx2[i],'members at h = '))[2]
      h=gsub(']','',h)
      lefNode=unlist(stringr::str_split(stx2[i],' members at h = '))[1]
      lefNode=unlist(stringr::str_split(lefNode,'branches and '))[2]
      sub1=paste0('C',i+1)
      all_node_obj[which(all_node_obj[,1]==paste0('C',i+1)),4]=paste0('C',i)
      #sub2=all_node_obj[which(all_node_obj[,4]==''),1]
      for(j in 1:nrow(all_node_obj)){
        #t.ind=which(all_node_obj[,1]==paste0('C',i+j+1))
        #if(length(t.ind))
        sub_mx=all_node_obj[which(all_node_obj[,1]==paste0('C',i+j+1)),]
        if(sub_mx[4]==''){
          sub2=sub_mx[1]
          break()
        }
      }
      all_node_obj[which(all_node_obj[,1]==sub2),4]=paste0('C',i)
      #all_node_obj=c(all_node_obj,list(Node=paste0('C',i),H=h,Members=lefNode,subNode1=paste0('C',i+1),subNode12=''))
      all_node_obj=rbind(all_node_obj,c(paste0('C',i),sub1,sub2,'',lefNode))
      nodeMap=rbind(nodeMap,c('',paste0('C',i),paste0('C',i+1),sub2,h))
    }else{
      node=gsub(' ','',unlist(stringr::str_split(stx2[i],'--leaf'))[2])
      node=gsub('^ \\s+','',node)
      node=gsub('\\s+ $','',node)
      node=gsub('^"','',node)
      node=gsub('"$','',node)
      #node=colnames(adjacency)[as.numeric(node)]
      leafs=rbind(leafs,c(node,paste0('C',i)))
      all_node_obj=rbind(all_node_obj,c(paste0('C',i),'','','',1))
      nodeMap=rbind(nodeMap,c(node,paste0('C',i),'','',0))
    }
  }
  #all_node_obj[all_node_obj[,4]=='',]
  #sort(leafs[,1])[1:10]
  od=rep(nrow(nodeMap),nrow(nodeMap))
  od[match(leafs[,1],nodeMap[,1])]=nrow(leafs):1
  nodeMap=nodeMap[order(od),]
  return(nodeMap)
}

getTreeNode_backup=function(hc_obj){
  dend=as.dendrogram(hc_obj)
  dend_label=partition_leaves(dend)
  dend_height=get_nodes_attr(dend,'height')
  nodeMap=rbind()
  for(i in length(dend_height):1){
    cls=dend_label[[i]]
    node1=paste0(cls,collapse = '\t')
    snode1=''
    snode2=''
    if(length(cls)>1){
      node1=openssl::md5(node1)
      snode1=paste0('C',i+1)
      node2_ts=dend_label[[i+1]]
      cls1=cls[which(!cls%in%node2_ts)]
      node2=paste0(cls1,collapse = '\t')
      if(length(cls1)>1){
        node2=openssl::md5(node2)
      }
      snode2=nodeMap[which(nodeMap[,1]==node2),2]
    }
    nodeMap=rbind(nodeMap,c(node1,paste0('C',i),snode1,snode2,dend_height[i]))
  }
  od=rep(10000000,nrow(nodeMap))
  od[match(hc_obj$labels[hc_obj$order],nodeMap[,1])]=1:length(hc_obj$labels)
  nodeMap=nodeMap[order(od),]
  nodeMap[which(nodeMap[,3]!=''),1]=''
  #head(nodeMap)
  #head(nodeMap1)
  #sum(nodeMap1[,1]==nodeMap1[,1])
  #nodeMap1[which(nodeMap1[,3]!=''),1]
  return(nodeMap)
}
getClusterTree=function(exp,dMethod='euclidean',hMethod='complete'){
  cNames=colnames(exp)
  colnames(exp)=paste0('C',1:ncol(exp))
  hc=hcluster(t(exp),method =dMethod,link=hMethod )
  nodeMap=getTreeNode(hc)
  nodeMap2=getTreeNode_backup(hc)
  #dim(nodeMap1)
  #head(nodeMap1)
  sum(nodeMap2[,2]==nodeMap[,2])
  sum(nodeMap2[,1]==nodeMap[,1])
  sum(nodeMap2[,4]==nodeMap[,4])
  cbind(nodeMap2[which(nodeMap2[,4]!=nodeMap[,4]),],nodeMap[which(nodeMap2[,4]!=nodeMap[,4]),])
  
  #dim(nodeMap)
  nodeMap[,1]=cNames[match(nodeMap[,1],colnames(exp))]
  nodeMap[is.na(nodeMap[,1]),1]=''
  #nodeMap[which(nodeMap[,3]!=''),1]
  #head(exp)
  #dend=as.dendrogram(hc)
  #dend_label=partition_leaves(dend)
  #dend_height=get_nodes_attr(dend,'height')
  #nodeMap=rbind()
  #for(i in length(dend_height):1){
  #  cls=dend_label[[i]]
  #  node1=paste0(cls,collapse = '\t')
  #  snode1=''
  #  snode2=''
  #  if(length(cls)>1){
  #    node1=openssl::md5(node1)
  #    snode1=paste0('C',i+1)
  #    node2_ts=dend_label[[i+1]]
  #    cls1=cls[which(!cls%in%node2_ts)]
  #    node2=paste0(cls1,collapse = '\t')
  #    if(length(cls1)>1){
  #      node2=openssl::md5(node2)
  #    }
  #    snode2=nodeMap[which(nodeMap[,1]==node2),2]
  #  }
  #  nodeMap=rbind(nodeMap,c(node1,paste0('C',i),snode1,snode2,dend_height[i]))
  #}
  #od=rep(10000000,nrow(nodeMap))
  #od[match(hc$labels[hc$order],nodeMap[,1])]=1:length(hc$labels)
  #nodeMap=nodeMap[order(od),]
  #nodeMap[which(nodeMap[,3]!=''),1]=''
  #head(nodeMap)
  #head(nodeMap1)
  #sum(nodeMap1[,1]==nodeMap1[,1])
  #nodeMap1[which(nodeMap1[,3]!=''),1]
  #head(nodeMap)
  return(nodeMap)
}

library(WGCNA)

if(method=='module'){

  exp_path=unlist(data$exp_path)
  geneMode=paste0(as.character(unlist(data$geneMode)),'')
  geneModeValue=unlist(data$geneModeValue)
  sampleCut=unlist(data$sampleCut)
  netType=unlist(data$netType)
  outFolder=opt$outfile
  minModuleSize=unlist(data$minModuleSize)
  deepSplit=unlist(data$deepSplit)
  oldFolder=unlist(data$oldFolder)
  mergeCutHeight=unlist(data$mergeCutHeight)
  
  geneModeValue=as.numeric(geneModeValue)
  
if(oldFolder!=''&oldFolder!=outFolder&file.exists(paste0(oldFolder,'/wgcna_TOM.RData'))){
  load(paste0(oldFolder,'/wgcna_TOM.RData'))
  file.copy(paste0(oldFolder,'/wgcna_pre_tree.txt'),paste0(outFolder,'/wgcna_pre_tree.txt'),overwrite = T)
  file.copy(paste0(oldFolder,'/wgcna_exp_powers.txt'),paste0(outFolder,'/wgcna_exp_powers.txt'),overwrite = T)
  file.copy(paste0(oldFolder,'/wgcna_exp_data.txt'),paste0(outFolder,'/wgcna_exp_data.txt'),overwrite = T)
  file.copy(paste0(oldFolder,'/wgcna_dissTOM_tree.txt'),paste0(outFolder,'/wgcna_dissTOM_tree.txt'),overwrite = T)
  file.copy(paste0(oldFolder,'/wgcna_module_net.fst'),paste0(outFolder,'/wgcna_module_net.fst'),overwrite = T)
  file.copy(paste0(oldFolder,'/wgcna_TOM.RData'),paste0(outFolder,'/wgcna_TOM.RData'),overwrite = T)
  file.copy(paste0(oldFolder,'/network_edge.txt'),paste0(outFolder,'/network_edge.txt'),overwrite = T)
  file.copy(paste0(oldFolder,'/network_node.txt'),paste0(outFolder,'/network_node.txt'),overwrite = T)
  logs=c(logs,paste0('SampleList=',paste0(row.names(filter.exp),collapse = '#,#')))
  logs=c(logs,paste0('saved module network data:',paste0(cytq,collapse = ',')))
  logs=c(logs,paste0('set power=',cutPower,',By sample count=',nrow(filter.exp)))
  logs=c(logs,paste0('saved TOM data'))
  #wgcna_module_net.fst
  
}else{

  
dat=data.table::fread(exp_path, sep = "\t",header = T,stringsAsFactors = F,check.names = F
                      ,na.strings="NA",data.table = F)  
dat=dat[match(unique(dat[,1]),dat[,1]),]
row.names(dat)=dat[,1]
dat=dat[,-1]
logs=c(logs,'filter SD=0')
sds=apply(dat, 1, sd)
dat=dat[which(sds>0),]
logs=c(logs,paste0('filtered SD=0,row=',nrow(dat),',col=',ncol(dat)))

if(geneMode=='mad'){
  logs=c(logs,paste0('filter MAD>',ceiling(geneModeValue),'%'))
  mads=apply(dat, 1, mad)
  t_inds=which(mads>=quantile(mads,seq(0,1,0.01))[ceiling(geneModeValue)+1])
  dat=dat[t_inds,]
  logs=c(logs,paste0('filtered MAD,row=',nrow(dat),',col=',ncol(dat)))
}else if(geneMode=='sd'){
  logs=c(logs,paste0('filter SD>',ceiling(geneModeValue),'%'))
  mads=apply(dat, 1, sd)
  t_inds=which(mads>=quantile(mads,seq(0,1,0.01))[ceiling(geneModeValue)+1])
  dat=dat[t_inds,]
  logs=c(logs,paste0('filtered SD,row=',nrow(dat),',col=',ncol(dat)))
}else if(geneMode=='mean'){
  logs=c(logs,paste0('filter mean>',(geneModeValue),'%'))
  mads=apply(dat, 1, mean)
  t_inds=which(mads>geneModeValue)
  dat=dat[t_inds,]
  logs=c(logs,paste0('filtered mean,row=',nrow(dat),',col=',ncol(dat)))
}
#logs=c()
#print('=====')
filter.exp=t(dat)
dim(filter.exp)
logs=c(logs,paste0('Sample clustering'))
sTree=getClusterTree(dat,hMethod = 'average')
write.table(sTree,file = paste0(outFolder,'/wgcna_pre_tree.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved sample clustered'))


gsg = goodSamplesGenes(filter.exp, verbose = 3);
#gsg$allOK
# if FALSE:
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    logs(logs,paste("Removing genes:", paste(colnames(filter.exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    logs(logs,paste("Removing samples:", paste(rownames(filter.exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  filter.exp = filter.exp[gsg$goodSamples, gsg$goodGenes]
}
logs=c(logs,paste0('SampleFiltered sampleNum=',nrow(filter.exp),'geneNum=',ncol(filter.exp)))

#if(sampleCut=='1'){
#  logs=c(logs,paste0('SampleFiltering'))
#  logs=c(logs,paste0('Sample clustering'))
#  sampleTree = hclust(dist(filter.exp), method = "average")
  
#  logs=c(logs,paste0('Sample clustered'))
#  q1=quantile(sampleTree$height, c(0,0.25, 0.5, 0.75,1), na.rm=T)
#  otl=1.5*(q1[4]-q1[2])
#  otl.u=ifelse(q1[4]+otl>q1[5],q1[5],q1[4]+otl)
#  otl.d=ifelse(q1[2]-otl<q1[1],q1[1],q1[2]-otl)
#  logs=c(logs,paste0('SampleFiltering,outlier upper=',otl.u))
  #outlier=which(sampleTree$height>=otl.d&sampleTree$height<=otl.u)
#  ctr=cutree(sampleTree,h = otl.u+1)
#  ctr.tb=table(ctr)
#  outlier=match(names(which(ctr==names(ctr.tb)[which.max(ctr.tb)])),row.names(filter.exp))
#  logs=c(logs,paste0('SampleFiltered sampleNum=',length(outlier)))
#  filter.exp=filter.exp[outlier,]
#}

write.table(cbind(Tag=colnames(filter.exp),t(filter.exp)),file = paste0(outFolder,'/wgcna_exp_data.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('SampleList=',paste0(row.names(filter.exp),collapse = '#,#')))

# pdf(file = paste0(outFolder,'/SampleCluster.pdf'),width = 12,height = 6)
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
# dev.off()
#height=140000
#print(height)
#print(min(sampleTree$height))

#print('===1')
logs=c(logs,paste0('run pickSoftThresholding'))
powers = c(c(1:10), seq(from = 12, to=30, by=2))
#print('===2')
#print(head(filter.exp))
sft = pickSoftThreshold((filter.exp), powerVector=powers, 
                        networkType=netType, verbose=5,blockSize = 20000)
cutPower=sft$powerEstimate
logs=c(logs,paste0('run pickSoftThresholded power=',cutPower))
if(is.na(cutPower)){
  if(netType=='unsigned'|netType=='signed_hybrid'){
    if(nrow(filter.exp)<20){
      cutPower=9
    }else if(nrow(filter.exp)<30){
      cutPower=8
    }else if(nrow(filter.exp)<40){
      cutPower=7
    }else{
      cutPower=6
    }
  }else{
    if(nrow(filter.exp)<20){
      cutPower=18
    }else if(nrow(filter.exp)<30){
      cutPower=16
    }else if(nrow(filter.exp)<40){
      cutPower=14
    }else{
      cutPower=12
    }
  }
}
logs=c(logs,paste0('set power=',cutPower,',By ',netType,' sample count=',nrow(filter.exp)))

sft.tab=sft$fitIndices
sft.tab$SFTR2=-sign(sft$fitIndices[,3])*sft$fitIndices[,2]
write.table(sft.tab,file = paste0(outFolder,'/wgcna_exp_powers.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved powerEstimateed'))


logs=c(logs,paste0('runing adjacency'))
adjacency = adjacency(filter.exp, power = cutPower,type = netType)
logs=c(logs,paste0('runing TOMsimilarity'))
TOM = TOMsimilarity(adjacency)
colnames(TOM)=colnames(adjacency)
row.names(TOM)=row.names(adjacency)
#head(TOM[,1:10])
dissTOM = 1-TOM
plotTOM = dissTOM^7
#dim(TOM)
diag(plotTOM) = NA

#dim(adjacency)
#colnames(plotTOM)
# Call the plot function
#colnames(dissTOM)=colnames(adjacency)
#row.names(dissTOM)=colnames(adjacency)
cName_TOM=colnames(dissTOM)
rName_TOM=row.names(dissTOM)
colnames(dissTOM)=paste0('C',1:ncol(dissTOM))
row.names(dissTOM)=paste0('R',1:nrow(dissTOM))
logs=c(logs,paste0('hclust dissTOM'))
geneTree = hclust(as.dist(dissTOM), method = "average")
dim(dissTOM)

pdf(file = paste0(outFolder, "/gene_clustering.pdf"))
plot(geneTree, main="Hierarchical Clustering of Genes", xlab="Genes", ylab="Distance")
dev.off()

#plot(geneTree)
logs=c(logs,paste0('hclusted dissTOM'))
nodeMap=getTreeNode(geneTree)
#head(nodeMap)
nodeMap[,1]=rName_TOM[match(nodeMap[,1],row.names(dissTOM))]
nodeMap[is.na(nodeMap[,1]),1]=''
colnames(dissTOM)=cName_TOM
row.names(dissTOM)=rName_TOM
#dend=as.dendrogram(geneTree)
#capture.output(str(dend))->stx2
#write.table(stx,row.names = F,col.names = F ,paste0(outFolder,'/wgcna_dissTOM_tree.text'))
#gsub('^ \\s+','',stx2)->stx2
#nodeMap=rbind()
#leafs=rbind()
#for(i in length(stx2):1){
#  i=length(stx2)
#  if(length(grep('--leaf',stx2[i]))==0){
#    h=unlist(stringr::str_split(stx2[i],'members at h = '))[2]
#    h=gsub(']','',h)
#    nodeMap=rbind(nodeMap,c('',paste0('C',i),paste0('C',i+1),paste0('C',i+2),h))
#  }else{
#    node=gsub(' ','',unlist(stringr::str_split(stx2[i],'--leaf'))[2])
#    node=gsub('^ \\s+','',node)
#    node=gsub('\\s+ $','',node)
#    node=gsub('^"','',node)
#    node=gsub('"$','',node)
    #node=colnames(adjacency)[as.numeric(node)]
#    leafs=rbind(leafs,c(node,paste0('C',i)))
#    nodeMap=rbind(nodeMap,c(node,paste0('C',i),'','',0))
#  }
#}
#sort(leafs[,1])[1:10]
#od=rep(nrow(nodeMap),nrow(nodeMap))
#od[match(leafs[,1],nodeMap[,1])]=nrow(leafs):1
#nodeMap=nodeMap[order(od),]
#head(nodeMap)
#geneTree$labels

#match(colnames(dissTOM)[geneTree$order],nodeMap[,1])

#which(nodeMap[,3]=='')

#print(geneTree$height)
#class(dend)
#summary(dend)
#attr(dend,x = length)
#plot(geneTree,hang = -1, cex = 0.6)
#plot(as.phylo(geneTree), type = "fan")
#dend_label=partition_leaves(dend)
#dend_height=get_nodes_attr(dend,'height')
#nodeMap=rbind()
#for(i in length(dend_height):1){
#  cls=dend_label[[i]]
#  node1=paste0(cls,collapse = '\t')
#  snode1=''
#  snode2=''
#  if(length(cls)>1){
#    node1=openssl::md5(node1)
#    snode1=paste0('C',i+1)
#    node2_ts=dend_label[[i+1]]
#    cls1=cls[which(!cls%in%node2_ts)]
#    node2=paste0(cls1,collapse = '\t')
#    if(length(cls1)>1){
#      node2=openssl::md5(node2)
#    }
#    snode2=nodeMap[which(nodeMap[,1]==node2),2]
#  }
#  nodeMap=rbind(nodeMap,c(node1,paste0('C',i),snode1,snode2,dend_height[i]))
#}
#nodeMap[which(nodeMap[,3]!=''),1]=''
#od=rep(10000,nrow(nodeMap))
#od[match(geneTree$labels[geneTree$order],nodeMap[,1])]=1:length(geneTree$labels)
#nodeMap=nodeMap[order(od),]
write.table(nodeMap,file = paste0(outFolder,'/wgcna_dissTOM_tree.txt'),quote = F,row.names = F,sep = '\t')
#head(nodeMap)
logs=c(logs,paste0('saved hclust dissTOM tree'))
logs=c(logs,paste0('exporting module network data'))
#colnames(TOM)
#row.names(TOM)
cyt = exportNetworkToCytoscape(TOM,threshold =0)
#head(cyt$edgeData)
tidyfst::export_fst(cyt$edgeData[,1:3],path = paste0(outFolder,'/wgcna_module_net.fst')) 
cytq=quantile(cyt$edgeData[,3])
logs=c(logs,paste0('saved module network data:',paste0(cytq,collapse = ',')))
save(filter.exp,adjacency,dissTOM,geneTree,cytq,cutPower,file = paste0(outFolder,'/wgcna_TOM.RData'))
logs=c(logs,paste0('saved TOM data'))
}

#load(paste0('/pub1/data/mg_projects/projects/web_script/tool_runing/97d971bc04475222ec72f4a05e66910f/wgcna_TOM.RData'))

#  plotDendroAndColors(geneTree,dynamicColors ,
#                      "Module colors",
#                      dendroLabels = FALSE, hang = 0.03,
#                      addGuide = TRUE, guideHang = 0.05)
#################cal module#################
#dim(dissTOM)
#tom_tree=readMatrix(paste0(oldFolder,'/wgcna_dissTOM_tree.txt'),row=F)
#head(tom_tree)
#setdiff(colnames(filter.exp),unique(tom_tree[,1]))[1:10]
#setdiff(unique(tom_tree[,1]),colnames(filter.exp))[1:10]

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM
                            ,deepSplit = as.numeric(deepSplit), pamRespectsDendro = FALSE,
                            minClusterSize = as.numeric(minModuleSize))
dynamicColors = labels2colors(dynamicMods)
#MEList = moduleEigengenes(filter.exp, colors = dynamicColors)
merge = mergeCloseModules(filter.exp, dynamicColors, cutHeight = as.numeric(mergeCutHeight), verbose = 3)
# The merged module colors
mergedColors = merge$colors
#dim(m_dat_sum)
#length(dynamicColors)
m_dat_sum=as.data.frame(table(mergedColors))
logs=c(logs,paste0('module count:',paste0(paste0(m_dat_sum[,1],'-',m_dat_sum[,2]),collapse = ',')))
write.table(m_dat_sum,file = paste0(outFolder,'/wgcna_module_summary.txt'),quote = F,row.names = F,sep = '\t')

#head(cbind(Gene=colnames(filter.exp),module=mergedColors))
write.table(cbind(Gene=colnames(filter.exp),dynamicModule=dynamicColors,mergeModule=mergedColors)
            ,file = paste0(outFolder,'/wgcna_module.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved hclust mergedColors data'))

MEs=merge$newMEs
colnames(MEs)=gsub('^ME','',colnames(MEs))
write.table(cbind(Tag=row.names(MEs),MEs),file = paste0(outFolder,'/wgcna_MEs.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved hclust ME data'))
MEDiss = 1-cor(MEs)

sTree2=getClusterTree(MEDiss,dMethod = 'pearson',hMethod = 'average')
write.table(sTree2,file = paste0(outFolder,'/wgcna_ME_tree.txt'),quote = F,row.names = F,sep = '\t')
m_inds=match(sTree2[which(sTree2[,1]!='')],colnames(MEDiss))
write.table(MEDiss[m_inds,m_inds],file = paste0(outFolder,'/wgcna_ME_tree_heat.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved hclust ME tree'))
#####################Module membership########################
geneModuleMembership = as.data.frame(cor(filter.exp, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(filter.exp)));
#melt(geneModuleMembership)
#colnames(geneModuleMembership)
#cbind(colnames(filter.exp),mergedColors)
mmr=reshape::melt(cbind(Tag=row.names(geneModuleMembership),geneModuleMembership),id.vars=c('Tag'))
mmp=reshape::melt(cbind(Tag=row.names(MMPvalue),MMPvalue),id.vars=c('Tag'))
mmr=mmr[paste0(mmr[,1],'\t',mmr[,2])%in%paste0(colnames(filter.exp),'\t',mergedColors),]
mmr$pvalue=mmp[match(paste0(mmr[,1],'\t',mmr[,2]),paste0(mmp[,1],'\t',mmp[,2])),3]
colnames(mmr)=c('Gene','Module','R','pvalue')
mmr=mmr[match(colnames(filter.exp),mmr[,1]),]
write.table(mmr,file = paste0(outFolder,'/wgcna_MM_value.txt'),quote = F,row.names = F,sep = '\t')
logs=c(logs,paste0('saved Module membership data'))
#mergedColors
if(!file.exists(paste0(outFolder,'/clini.json'))){
  write.table('',file = paste0(outFolder,'/clini.json'),quote = F,row.names = F,sep = '\t')
}
tidyfst::export_fst(as.data.frame(cbind(Gene=colnames(filter.exp),mergeModule=mergedColors)),path = paste0(outFolder,'/wgcna_module.fst')) 
tidyfst::export_fst(as.data.frame(cbind(Sample=row.names(MEs),MEs)),path = paste0(outFolder,'/wgcna_MEs.fst')) 
tidyfst::export_fst(mmr,path = paste0(outFolder,'/wgcna_MM_value.fst')) 
tidyfst::export_fst(as.data.frame(cbind(Sample=row.names(filter.exp),filter.exp)),path = paste0(outFolder,'/filter.exp.fst')) 
logs=c(logs,paste0('saved all module data'))

}else if(method=='cli'){
  
  outFolder=opt$outfile
  cNames=unlist(data$cNames)
  cTypes=unlist(data$cTypes)
  cliDatas=(data$cliDatas)
  oldFolder=unlist(data$oldFolder)
  #cliDatas[1][[2]]
  
  cliDatas=cliDatas[[1]]
  
  #table(cliDatas[,4])
  
  file.copy(opt$infile,paste0(oldFolder,'/clini.json'),overwrite = T)
  #data<-jsonlite::stream_in(file(paste0(oldFolder,'/clini.json')),pagesize = 1000)
  
  #oldFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/97d971bc04475222ec72f4a05e66910f'
  
  baseFolder=oldFolder
  #geneModule=tidyfst::import_fst(paste0(baseFolder,'/wgcna_module.fst'),as.data.table = F)
  mes=tidyfst::import_fst(paste0(baseFolder,'/wgcna_MEs.fst'),as.data.table = F)
  #head(mes)
  row.names(mes)=mes[,1]
  mes=mes[,-1]
  cmp=intersect(row.names(mes),cliDatas[,1])
  mes=mes[match(cmp,row.names(mes)),]
  cliDatas=cliDatas[match(cmp,cliDatas[,1]),]
  nColNames=c()
  nColSubNames=c()
  nColData=cbind()
  for(i in 1:length(cTypes)){
    if(cTypes[i]>0){
      nColNames=c(nColNames,cNames[i])  
      nColSubNames=c(nColSubNames,cNames[i])
      nColData=cbind(nColData,cliDatas[,i+1])
    }else{
      dt=cliDatas[,i+1]
      udt=unique(dt)
      udt=udt[which(!is.na(udt)&udt!='')]
      for(u in udt){
        ndt=ifelse(dt==u,1,0)
        ndt[which(is.na(dt))]=NA
        nColData=cbind(nColData,ndt)
        nColSubNames=c(nColSubNames,u)
        nColNames=c(nColNames,cNames[i])
      }
    }
  }
  nColData=apply(nColData, 2, as.numeric)
  #nColData=nColData[,1]
  cliModuleMembership = as.data.frame(cor(mes,nColData, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(cliModuleMembership), nrow(mes)));
  colnames(cliModuleMembership)=paste0(nColNames,'(',nColSubNames,')')
  colnames(MMPvalue)=paste0(nColNames,'(',nColSubNames,')')
  
  cNs=rbind(c(-1,'',nColNames),c(-1,'',nColSubNames))
  colnames(cNs)=c('isPvalue','Module',colnames(cliModuleMembership))
  cdt2=rbind(cbind(isPvalue=0,Module=row.names(cliModuleMembership),cliModuleMembership),
        cbind(isPvalue=1,Module=row.names(MMPvalue),MMPvalue))
  cdt2=apply(cdt2, 2,as.character)
  
  write.table(rbind(as.data.frame(cNs),as.data.frame(cdt2)),file = paste0(baseFolder,'/wgcna_ME_cli_cor.txt')
              ,quote = F,row.names = F,sep = '\t')
  
  #rbind(cliModuleMembership,MMPvalue)
  #geneModule=tidyfst::import_fst(paste0(baseFolder,'/wgcna_module.fst'),as.data.table = F)
  mmr=tidyfst::import_fst(paste0(baseFolder,'/wgcna_MM_value.fst'),as.data.table = F)
  filter.exp=tidyfst::import_fst(paste0(baseFolder,'/filter.exp.fst'),as.data.table = F)
  row.names(filter.exp)=filter.exp[,1]
  filter.exp=filter.exp[,-1]
  filter.exp=filter.exp[match(cmp,row.names(filter.exp)),]
  #head(filter.exp[,1:10])
  
  cliGeneMembership = as.data.frame(cor(filter.exp,nColData, use = "p"));
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(cliGeneMembership), nrow(mes)));
  colnames(cliGeneMembership)=paste0(nColNames,'(',nColSubNames,')')
  colnames(MMPvalue)=paste0(nColNames,'(',nColSubNames,')')
  #head(gs)
  gs=cbind(mmr,cliGeneMembership[match(mmr[,1],row.names(cliGeneMembership)),]
           ,MMPvalue[match(mmr[,1],row.names(MMPvalue)),])
  write.table(gs,file = paste0(baseFolder,'/wgcna_MM_GS_cor.txt')
              ,quote = F,row.names = F,sep = '\t')  
  if(ncol(cliGeneMembership)==1){
    ngm=cbind(abs(cliGeneMembership[match(mmr[,1],row.names(cliGeneMembership)),]))
    colnames(ngm)=colnames(cliGeneMembership)
  }else{
    ngm=abs(cliGeneMembership[match(mmr[,1],row.names(cliGeneMembership)),])
  }
  #head(ngm)
  all_cr=rbind()
  for(u in unique(as.character(mmr[,2]))){
    #u=unique(as.character(mmr[,2]))[1]
    t_inds=which(mmr[,2]==u)
    if(length(t_inds)>3){
      cr=as.data.frame(cor(cbind(abs(mmr[t_inds,3])),ngm[t_inds,], use = "p"));
      cr_p=as.data.frame(corPvalueStudent(as.matrix(cr), length(t_inds)));
      cdt=cbind(M=u,Clinical=1:ncol(ngm),t(rbind(cr,cr_p)))
      colnames(cdt)=c('Module','Clinical','R','P')
      all_cr=rbind(all_cr,cdt)
    }else{
      cdt=cbind(rep(u,ncol(ngm)),1:ncol(ngm),rep(NA,ncol(ngm)),rep(NA,ncol(ngm)))
      colnames(cdt)=c('Module','Clinical','R','P')
      all_cr=rbind(all_cr,cdt)
    }
  }
  write.table(all_cr,file = paste0(baseFolder,'/wgcna_MM_GS_cor_RP.txt'),quote = F,row.names = F,sep = '\t') 
  
}else if(method=='export'){
  
  weight=as.numeric(unlist(data$weight))
  modules=unlist(data$modules)
  fileTree=as.numeric(unlist(data$fileTree))
  baseFolder=unlist(data$oldFolder)
  print(weight)
  print(fileTree)
  #cliDatas[1][[2]]
  #table(cliDatas[,4])
  
  #data<-jsonlite::stream_in(file(paste0(oldFolder,'/export.json')),pagesize = 1000)
  
  #oldFolder='/pub1/data/mg_projects/projects/web_script/tool_runing/97d971bc04475222ec72f4a05e66910f'
  if(fileTree==1){
    al.dt=c()
    for(fl in dir(oldFolder)){
      if(length(grep('.json$',fl))>0|length(grep('.log$',fl))>0){
      }else{
        al.dt=c(al.dt,fl)    
      }
    }
    write.table(al.dt,file = paste0(opt$outfile,'/file.tree'),quote = F,row.names = F,col.names = F,sep = '\t')
  }else{
    file.copy(opt$infile,paste0(oldFolder,'/export.json'),overwrite = T)
    #weight=0.2
    #modules=c('darkgreen','turquoise');
    net=tidyfst::import_fst(paste0(baseFolder,'/wgcna_module_net.fst'),as.data.table = F)
    net=net[net[,3]>weight,]
    #dim(net)
    geneModule=tidyfst::import_fst(paste0(baseFolder,'/wgcna_module.fst'),as.data.table = F)
    genes=as.character(geneModule[which(as.character(geneModule[,2])%in%modules),1])
    t1.inds=which(as.character(net[,1])%in%genes)
    t2.inds=which(as.character(net[,2])%in%genes)
    t.inds=intersect(t1.inds,t2.inds)
    if(length(t.inds)==0){
      write.table(rbind(colnames(net)),file = paste0(baseFolder,'/network_edge.txt'),quote = F,row.names = F
                  ,sep = '\t',col.names = F)
      write.table(rbind(colnames(geneModule)),file = paste0(baseFolder,'/network_node.txt'),quote = F,row.names = F
                  ,sep = '\t',col.names = F)
    }else{
      net=net[t.inds,]
      ug=unique(c(as.character(net[,1]),as.character(net[,2])))
      node=geneModule[match(ug,as.character(geneModule[,1])),]
    
      write.table(net,file = paste0(baseFolder,'/network_edge.txt'),quote = F,row.names = F,col.names = T,sep = '\t')
      write.table(node,file = paste0(baseFolder,'/network_node.txt'),quote = F,row.names = F,col.names = T,sep = '\t')
    }
  }
  
}

},error = function(e) {
  print(conditionMessage(e))
  logs=c(logs,paste0('error:',conditionMessage(e)))
}, finally = {
  write.table(logs,file = paste0(opt$outfile,'/run.log'),quote = F,row.names = T,col.names = T,sep = '\t')
})


