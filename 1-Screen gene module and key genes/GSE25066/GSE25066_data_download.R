library(BiocManager)
#BiocManager::install("GEOquery")
#BiocManager::install("DOSE")
library(GEOquery)
library(DOSE)
setwd("D:\\R_WorkSpace\\GSE25066")
options( 'download.file.method.GEOquery' = 'libcurl' )
eSet<-getGEO("GSE25066",destdir=".",#下载在当前目录
             getGPL=T)#平台信息不要
#提取表达矩阵
exprSet <-exprs(eSet[[1]])  #此时，基因表达矩阵是探针和样本

fdata = fData(eSet[[1]])#需要getGPL = T,平台信息,注释文件包含探针与基因的对应关系
pdata = pData(eSet[[1]]) #表型信息，里面包含了需要的分组信息
samplenames <- sampleNames(eSet[[1]])
dim(exprSet )
dim(fdata)
dim(pdata)
length(samplenames)
#找出ID对应的基因
genes <- fdata[,c("ID","Gene Symbol")]
exprSet <- data.frame(exprSet )
exprSet$ID <- rownames(exprSet)
exprSet <- merge(exprSet,genes,by.x = "ID",by.y = 'ID')
exprSet <- exprSet[,-1]
dim(exprSet)
#处理重复基因
rowMeans = apply(exprSet[,c(1:12)],1,function(x) mean(as.numeric(x), na.rm = T))####计算行平均值
rowMeans_2 = data.frame(rowMeans)
exprSet = exprSet[order(rowMeans, decreasing = T),] #排序
exprSet_2 = exprSet[!duplicated(exprSet[, dim(exprSet)[2]]),] #去掉重复基因
exprSet_na = na.omit(exprSet_2)   #去掉缺失值
explan_final = exprSet_na[exprSet_na$`Gene Symbol`!= "",]
explan_final <- data.frame(explan_final[-grep("/",explan_final$`Gene Symbol`),]) #去一对多，grep是包含的意思，-就是不包含           

rownames(explan_final) <- explan_final$Gene.Symbol
explan_final<-explan_final[,-ncol(explan_final)]
write.csv(explan_final,'GSE25066_Data.csv',row.names = TRUE)
