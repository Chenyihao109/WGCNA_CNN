#install.packages("rjson")
#install.packages("jsonlite")
#整理tcga测试数据
setwd("D:\\R_WorkSpace\\Gene_regulatinh networks\\Breast_data")  #数据路径

library("rjson")
json <- jsonlite::fromJSON("metadata.cart.2022-12-28.json")   #metadata文件名，这个文件存了每一个文件对应的样本名
View(json)

sample_id <- sapply(json$associated_entities,function(x){x[,1]})#取出json文件下面提取出json文件的associated_entities列中的第一个元素
#sample_id就是所有的样本名
file_sample <- data.frame(sample_id,file_name=json$file_name)#得到文件名和样本名的对应 

#获取gdc_download文件夹下的所有TSV表达文件的 路径+文件名，这里路径要改为自己的路径
count_file <- list.files('gdc_download_20221228_062657.757930',pattern = '*.tsv',recursive = TRUE)  #Counts文件夹名"
count_file
#在count_file中分割出文件名
count_file_name <- strsplit(count_file,split='/')
#count_file_name
count_file_name <- sapply(count_file_name,function(x){x[2]})
count_file_name   #记录所有的文件名


#下面的修改基因数
matrix = data.frame(matrix(nrow=60660,ncol=0))

#下面的修改样本例数
for (i in 1:length(count_file_name)){
  path = paste0("gdc_download_20221228_062657.757930/",count_file[i])  #count文件夹名
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  
  colnames(data)<-data[2,]
  
  data <-data[-c(1:6),]
  
  data <- data[7]
  #数据类型,选择其中之一 3:unstranded;4:stranded_first；5：stranded_second；6：tpm_unstranded；7：fpkm_unstranded；8：fpkm_uq_unstranded
  
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]
  
  matrix <- cbind(matrix,data)
}
#matrix就是count矩阵，行名是ensID,列名是样本，这里包含了所有类型的基因，下一步筛选protein_coding
write.csv(matrix,'GBM FPKM Counts_matrix.csv',row.names = TRUE)#如果使用fpk，输出文件改为GBM FPKM_matrix.csv


##------------------------------增加部分：设置Gene Symbol为列名的矩阵（前面得到的matrix是Ensembl ID）------------------------------------------
#将ens转变为symbol
path=paste0("gdc_download_20221228_062657.757930/",count_file[1])
count_file[1]
data<-as.matrix(read.delim(path,fill=TRUE,header = FALSE,row.names = 1))
View(data)

#gene_name
gene_name<-data[-c(1:6),1]
#gene_type
gene_type<-data[-c(1:6),2]
table(gene_type)

matrix0<-cbind(gene_name,gene_type,matrix) #将基因名,基因类型和表达矩阵拼接
View(matrix0)

#筛选protein_coding
mRNA_matrix <- matrix0[matrix0$gene_type == "protein_coding",] ; 
dim(mRNA_matrix) # 拆分得到mRNA的表达矩阵 mRNA_matrix
View(mRNA_matrix)
#将gene_name列去除重复的基因，保留每个基因最大表达量结果

matrix0<-aggregate(.~gene_name,data=mRNA_matrix,max)


#将gene_name列设为行名
rownames(matrix0)<-matrix0[,1]
matrix0<-matrix0[,-c(1,2)]   #去掉冗余列
View(matrix0)
#write.csv(matrix0,"FPKM_count_matrix_name.csv",row.names = TRUE)
write.csv(matrix0,"protein_coding_gene_expression_matrix.csv",row.names = TRUE)



#************************************临床数据整理********************************
entity_submitter_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
sample_case <- t(rbind(entity_submitter_id,case_id))

clinical <- read.delim('clinical.cart.2022-12-28\\clinical.tsv',header = T)
clinical <- as.data.frame(clinical[duplicated(clinical$case_id),])

clinical_matrix <- merge(sample_case,clinical,by="case_id",all.x=T)
clinical_matrix <- clinical_matrix[,-1]
clinical_matrix2 <- clinical_matrix[,c("entity_submitter_id","vital_status","days_to_last_follow_up")]

write.csv(clinical_matrix,'clinical_matrix.csv',row.names = TRUE)
