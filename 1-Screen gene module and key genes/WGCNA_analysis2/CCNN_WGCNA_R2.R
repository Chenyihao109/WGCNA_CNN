
#options()$repos
#options()$BioC_mirror
#options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#options("repos"=c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
#options()$repos
#options()$BioC_mirror

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install('WGCNA',force = TRUE)
#BiocManager::install('org.Mm.eg.db')
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("org.Mm.eg.db")
#
library(org.Mm.eg.db)
options(stringsAsFactors = F)   #读取数据时，不将字符型转为因子
Sys.setenv(LANGUAGE="en")
library(BiocManager)
library(FactoMineR)

#BiocManager::install('factoextra')
library(factoextra)    
#BiocManager::install('tidyverse')
library(tidyverse)    #ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
library(data.table)   #多核读取文件
#installed.packages()[, c('Package', 'LibPath')]

#BiocManager::install('AnnotationDbi')

#BiocManager::install("Biobase")
#BiocManager::install('Biostrings')
#BiocManager::install('GenomeInfoDbData')
#BiocManager::install('GO.db')
library("GO.db")
#BiocManager::valid("GO.db")
library(WGCNA)
### 启用WGCNA多核计算
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores())
corType="pearsoon"
setwd("D:\\R_WorkSpace\\Gene_regulatinh networks\\WGCNA_analysis2")
#获取表达矩阵   行是基因名  列是样本名
exp=read.csv("protein_coding_gene_expression_matrix.csv",check.names = FALSE)
#view(exp)
rownames(exp)<-exp[,1]    #将行名改为基因名
exp=exp[,-1]              #删除多余的基因列
exp[1,]
#判断每个样本是肿瘤还是正常样本
sample <- colnames(exp)  #样本名
sample[1:4]
trait <- c()
trait_NO <- c()
for (i in 1:length(sample)){
  if(as.integer(substring(sample[i],14,15))>10){    #14、15位置大于10的为normal样本
    trait <- append(trait,"normal")
    trait_NO <- append(trait_NO,1)
  } else {
    trait <- append(trait,"tumor")
    trait_NO <- append(trait_NO,0)
  }
}
table(trait)
trait_NO[1:4]
datTraits <- data.frame(row.names = colnames(exp),group=trait,group_NO=trait_NO)
#### 筛选MAD前12500的基因
exp1<-log2(exp+1)    #取对数是为了减轻表达差异
keep_data<-exp1[order(apply(exp1,1,mad),decreasing = T)[1:12500],]   #计算每个基因在所有样本表达值的mad值
### 转置
datExpr <- as.data.frame(t(keep_data))   #筛选的10000个基因的表达矩阵（行基因 列样本）转置为（行样本，列基因）

##选择关键基因差异分析的对应矩阵
differenceAnalysis_matrix=datExpr
differenceAnalysis_matrix=cbind(trait,differenceAnalysis_matrix)
dffferenceAnalysis_matrix<- differenceAnalysis_matrix[,grepl("trait|FOXM1|E2F1|CCNB1|UBE2C|CKS2|NEK2|CDK1|CDCA8|TPX2|CEP55", colnames(differenceAnalysis_matrix))] #提取列
write.csv(dffferenceAnalysis_matrix,'differenceAnalysis_matrix.csv',row.names = TRUE)

save(datExpr,datTraits,file = "exp_and_group.Rdata")

############################## 1.判断数据质量 ################################
rm(list=ls())
load("exp_and_group.Rdata")
### 判断数据质量--缺失值
gsg <- goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
gsg$allOK
##################################################.聚类做离群样本检测
#我们是对离群样本做检测，目标是样本，所以行标得是样本名，故用我们上面处理得到矩阵datExpr0
#目的是通过样本的聚类树看看是否有任何明显的异常值
#dist(datExpr)            #dist()就是求样本与样本之间的距离，因为datExpr的行是样本
if(T){
  #针对样本做聚类树
  sampleTree <- hclust(dist(datExpr), method = "average")  #用距离作为样本之前的相关性，进行聚类
  par(mar = c(0,5,2,0))## 设置图形的边界，下，左，上，右的页边距 #par(cex = 0.6)	# 控制图片中文字和点的大小
  #画树形图
  plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 2,
       cex.axis = 1, cex.main = 1,cex.lab=1)
  ## 若样本有性状、表型，可以添加对应颜色，查看是否聚类合理
  #group记录样本的分组，这里有问题，下次做的时候，记得分类
  sample_colors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                                  colors = rainbow(length(table(datTraits$group))), 
                                  signed = FALSE)
  ## 绘制样品的系统聚类树及对应性状
  par(mar = c(1,4,3,1),cex=0.8)   #mar就是图片布局，四边的留白，cex设置标记大小
  pdf("step1_Sample dendrogram and trait.pdf",width = 8,height = 6) #创建一个Pdf文档
  plotDendroAndColors(sampleTree, sample_colors,
                      groupLabels = "trait",
                      cex.dendroLabels = 0.8,
                      marAll = c(1, 4, 3, 1),
                      cex.rowText = 0.01,
                      repel = TRUE, #标签不重叠
                      main = "Sample dendrogram and trait" )  #在pdf上画图
  ## Plot a line to show the cut
  # abline(h = 23500, col = "red") #根据实际情况而定
  dev.off()  #关闭当前绘图设备
}
##若存在显著离群点；剔除掉
if(T){
  clust <- cutreeStatic(sampleTree, cutHeight = 135, minSize = 10) # cutHeight根据实际情况而定
  ## 以高度23500切割，要求留下的最少为10个
  table(clust)## 查看切割后形成的集合
  keepSamples <- (clust==1)
  datExpr <- datExpr[keepSamples, ]    # 将clust序号为1的放入keepSamples
  datTraits <- datTraits[keepSamples,]  ## ncol，crow分别表示提取矩阵的列数和行数
  dim(datExpr) 
}
save(datExpr,datTraits,file="exp_and_grouplist_2.Rdata")  #剔除了离群值
##将处理后的表达矩阵输出到文件中，后续使用python进行分析
#write.csv(datExpr,'final_expression_matrix.csv',row.names = TRUE)   #这一步最后我只保存筛选出的基因所对应的表达矩阵，并不是所有基因
write.csv(datTraits,'final_samle_trait.csv',row.names = TRUE)
### 判断数据质量 : PCA进行分组查看，就是利用Pca去查看各个样本的分组情况，是否分开
rm(list = ls())  
load("exp_and_grouplist_2.Rdata")   #记录表达矩阵（转置后的）还有样本的分组情况
group_list <- datTraits$group
dat.pca <- PCA(datExpr, graph = F) #datExpr行是样本，列是基因
pca <- fviz_pca_ind(dat.pca,
                    title = "Principal Component Analysis",
                    legend.title = "Groups",
                    geom.ind = c("point"), #"point","text"  这里只选Point只显示点，不显示名字
                    pointsize = 2,
                    labelsize = 4,
                    repel = TRUE, #标签不重叠
                    col.ind = group_list, # 分组上色
                    addEllipses = TRUE, # 添加置信椭圆
                    axes.linetype=NA,  # remove axeslines
                    mean.point=F#去除分组中心点
) +
  theme(legend.position = "none")+  # "none" REMOVE legend
  coord_fixed(ratio = 1) #坐标轴的纵横比
pca
ggsave(pca,filename= "step1_Sample PCA analysis.pdf", width = 8, height = 8)  #将结果保存到Pdf中

##保存数据
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
write.csv(datExpr, "dataProcessing_expression_matrix.csv")
save(nGenes,nSamples,datExpr,datTraits,file="pca_inspecting_after.Rdata")  #记录gene个数，样本个数，转置后的表达矩阵，样本分组情况


################################ 2.挑选最佳阈值power ###################################
#挑选阈值还是需要表达矩阵的转置，也就是样本作为行名
rm(list=ls())
load("exp_and_grouplist_2.Rdata")
R.sq_cutoff=0.8  #设置R^2 cut-off,默认为0.85
if(T){
  #设置power参数选择范围
  powers<-c(seq(1,20,by=1),seq(22,30,by=2))
  sft<-pickSoftThreshold(datExpr,networkType = "unsigned",
                         powerVector = powers, #这里就是软阈值向量
                         RsquaredCut = R.sq_cutoff,  #设置R^2 cut-off（默认为0.85）
                         verbose = 5)
  sft[[1]]#就是在软阈值集合中对应的最佳的取值
  sft[[2]] #记录每一个软阈值对应的情况
  #a=sft$fitIndices   #aft$fitIndices就是dataframe形式的sft[[2]]
  pdf("step2_power-value.pdf",width = 16,height = 12)
  # Plot the results: 寻找拐点，确认最终power取值
  par(mfrow = c(1,2));#mfrow = c(1,2)设置一行两列的布局
  cex1 = 1.1;
  # 作为软阈值功率函数的无标度拓扑拟合指数
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",cex.axis=1.3,cex.lab=1.3)
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red",cex.axis=2,cex.lab=2)
  # 这条线对应于使用h的R^2截止值
  abline(h=R.sq_cutoff ,col="red")
  
  # 作为软阈值功率函数的平均连通性
  plot(sft$fitIndices[,1], sft$fitIndices[,5],     #1列就是不同的软阈值，第5列就是mean
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",cex.axis=1.3,cex.lab=1.3)
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
  #R语言中text()函数同abline()函数，lines()函数一样属于低水平函数，即在已有绘图中添加相关图形。
  #text()函数的作用是在给定的x和y坐标的位置添加字符串。
  abline(h=100,col="red")
  dev.off()
}
sft$powerEstimate  #查看估计的最佳power
power = sft$powerEstimate

# 若无向网络在power小于15或有向网络power小于30内，没有一个power值使
# 无标度网络图谱结构R^2达到0.8且平均连接度在100以下，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if(is.na(power)){
  # 官方推荐 "signed" 或 "signed hybrid"
  # 为与原文档一致，故未修改
  type = "unsigned"
  nSamples=nrow(datExpr)
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))      
                 )
  )
}

save(sft, power, file = "power_value.Rdata")
##################### 3.一步法构建加权共表达网络，识别基因模块 ####################
#一步法构建加权共表达网络，识别基因模块
#主要调整参数为minModuleSize和mergeCutHeight , 数值越大模块越少
rm(list=ls())
load(file="exp_and_grouplist_2.Rdata")
load(file="power_value.Rdata")
if(T){
  net<-blockwiseModules(
    datExpr,
    power = power,
    maxBlockSize = ncol(datExpr),
    corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
    networkType = "unsigned",
    TOMType = "unsigned", 
    minModuleSize = 30,    ##越大模块越少
    mergeCutHeight = 0.25, ##越大模块越少
    numericLabels = TRUE, 
    saveTOMs = F,
    verbose = 3
  )
  table(net$colors)
  # power: 上一步计算的软阈值
  # maxBlockSize:计算机能处理的最大模块的基因数量(默认5000),16G内存可以处理2万个，
  # 计算资源允许的情况下最好放在一个block里面。
  # corType：计算相关性的方法；可选pearson(默认)，bicor。后者更能考虑离群点的影响。
  # networkType:计算邻接矩阵时，是否考虑正负相关性；默认为"unsigned",可选"signed", "signed hybrid"
  # TOMType：计算TOM矩阵时，是否考虑正负相关性；默认为"signed",可选"unsigned"。但是根据幂律转换的邻接矩阵(权重)的非负性，所以认为这里选择"signed"也没有太多的意义。
  # numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
  # saveTOMs：最耗费时间的计算，可存储起来供后续使用，
  # mergeCutHeight: 合并模块的阈值，越大模块越少,一般为0.25
  # minModuleSize: 每个模块里最少放多少个基因，设定越大模块越少
  # 输出结果根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
  # **0 (grey)**表示**未**分入任何模块的基因。
}
### 模块可视化，层级聚类树展示各个模块
if(T){
  #将标签转换为用于打印的颜色
  moduleColors<-labels2colors(net$colors)
  table(moduleColors)
  ##绘制树状图和下面的模块颜色
  pdf("step3_genes-modules_ClusterDendrogram.pdf",width = 16,height = 12)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  dev.off()
  save(net, moduleColors, file = "genes_modules.Rdata")
}

###########绘制共表达网络##################################
rm(list=ls())
load(file="exp_and_grouplist_2.Rdata")
load(file="power_value.Rdata")
load(file="genes_modules.Rdata")
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
#计算TOM矩阵
#dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
#save(dissTOM, file = "dissTOM.Rdata")
load(file="dissTOM.Rdata")
#一般不用所有基因画图，因为太大了
#plotTOM = dissTOM^7;
#diag(plotTOM) = NA;
#sizeGrWindow(9,9)
#TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#选400个基因画图
nSelect = 400
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
sizeGrWindow(9,9)
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
#改变热图的深色背景为白色背景：
library(gplots) # 需要先安装这个包才能加载。
myheatcol = colorpanel(250,'red',"orange",'lemonchiffon')
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=myheatcol)


###############################选模块，判断乳腺癌关键基因在哪个模块##############
rm(list=ls())
load(file = "genes_modules.Rdata")
load(file = "pca_inspecting_after.Rdata")
GeneToColor <- data.frame(row.names = colnames(datExpr),gene=colnames(datExpr),color=moduleColors)
GeneToColor["BRCA1",2]
GeneToColor["FOXM1",2]
GeneToColor["E2F1",2]
GeneToColor["UBE2C",2]
GeneToColor["TPX2",2]
GeneToColor["CEP55",2]
GeneToColor["NEK2",2]
GeneToColor["CDCA8",2]
GeneToColor["CDK1",2]
GeneToColor["CCNB1",2]
GeneToColor["CKS2",2]


yellowmodule<-GeneToColor[GeneToColor$color=="yellow",]
blackmodule<-GeneToColor[GeneToColor$color=="black",]
write.csv(yellowmodule, "yellowmodule.csv")
write.csv(blackmodule, "blackmodule.csv")

#################对感兴趣的模块进行GO富集分析####################################
rm(list=ls())
load(file="exp_and_grouplist_2.Rdata")
load(file="power_value.Rdata")
load(file="genes_modules.Rdata")
### 条件设置
OrgDb = "org.Hs.eg.db"  # "org.Mm.eg.db"  "org.Hs.eg.db"
genetype = "SYMBOL"    # "SYMBOL"   "ENSEMBL"
table(moduleColors)
choose_module <- c("black","yellow","brown")
#choose_module <- c("yellow")

if(T){
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  
  gene_module <- data.frame(gene=colnames(datExpr),
                            module=moduleColors)
  write.csv(gene_module,file = "step6_gene_moduleColors.csv",row.names = F, quote = F) 
  tmp <- bitr(gene_module$gene,fromType = genetype,  # "SYMBOL"   "ENSEMBL"
              toType = "ENTREZID",
              OrgDb = "org.Hs.eg.db")
  gene_module_entrz <- merge(tmp,gene_module, by.x=genetype, by.y="gene")
  
  choose_gene_module_entrz <- gene_module_entrz[gene_module_entrz$module %in% choose_module,]
  
  ###run go analysis
  formula_res <- compareCluster(
    ENTREZID~module,
    data = choose_gene_module_entrz,
    fun = "enrichGO",
    OrgDb = "org.Hs.eg.db",
    ont = "BP",  #One of "BP", "MF", and "CC"  or "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.25
  )
  
  ###精简GO富集的结果,去冗余
  lineage1_ego <- simplify( 
    formula_res,
    cutoff=0.5,
    by="p.adjust",
    select_fun=min
  )
  save(gene_module, formula_res, lineage1_ego, file="step6_module_GO_term.Rdata")
  write.csv(lineage1_ego@compareClusterResult,
            file="step6_module_GO_term.csv")
  
  ### 绘制dotplot图
  dotp <- dotplot(lineage1_ego,
                  showCategory=10,
                  includeAll = TRUE, #将有overlap的结果也展示出来
                  label_format=90)
  dotp
  ggsave(dotp,filename= "step6_module_GO_term.pdf", #device = cairo_pdf,
         width = 12, 
         height = 15)
}

####################### 4.关联基因模块与表型 #####################################
#绘制模块与表型的相关性热图、模块与表型的相关性boxplot图、基因与模块、表型的相关性散点图； 
#结合所得各类相关性结果，判断出哪些模块与表型是高度相关的，从而关联基因模块与表型。
#这里示范在绘制基因与模块、表型的相关性散点图时，选择要查看的离散性状的表型是“primed”
rm(list = ls())  
load(file = "pca_inspecting_after.Rdata")
load(file = "power_value.Rdata")
load(file = "genes_modules.Rdata")
length(table(moduleColors))
## 模块与表型的相关性热图
if(T){
  datTraits$group <- as.factor(datTraits$group)   #将每个样本的分组转为因子
  design <- model.matrix(~0+datTraits$group)  #构建分组矩阵 也就是差异分析的设计矩阵，将样本的表型转为one-hot编码
  colnames(design) <- levels(datTraits$group) #设置分组矩阵的列名
  MES0 <- moduleEigengenes(datExpr,moduleColors)$eigengenes  #计算模块的特征基因，这里使用的是一步到位的共表达网络的结果
  MEs <- orderMEs(MES0)  #将相近的特征向量放在一起
  moduleTraitCor <- cor(MEs,design,use = "p") #计算模块和表型的相似性
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  #计算pvalue值
  textMatrix <- paste0(signif(moduleTraitCor,2),"\n(",
                       signif(moduleTraitPvalue,1),")")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  pdf("step4_Module-trait-relationship_heatmap.pdf",  #画图
      width = 2*length(colnames(design)), 
      height = 0.6*length(names(MEs)) )
  par(mar=c(5, 12, 3, 3)) #留白：下、左、上、右
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels=c("tumor","normal"),
                 #xLabels = colnames(design),
                 yLabels = names(MEs),
                 ySymbols = names(MEs),
                 colorLabels = F,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = F,
                 cex.text = 0.7,
                 zlim = c(-1,1), 
                 xColorWidth = 0.5,
    
                 main = "Module-trait relationships")
  dev.off()
  save(design, file = "step4_design.Rdata")
}

table(moduleColors)

######将与肿瘤显著相关的模块所包含的基因筛选出来，构建蛋白质交互网络，作为CNN的先验网络
gene_name=colnames(datExpr)
length(gene_name)
choose_gene=c()
for(i in 1:ncol(datExpr)){
  if(moduleColors[i] %in% c("yellow")){
    choose_gene=append(choose_gene,gene_name[i])
  }
}
length(choose_gene)
table(moduleColors)
write.table(choose_gene,'yellowModule_primal_gene_list.txt',row.names = F,col.names = F,quote = F) #字符不要带引号
datExpr2=datExpr[,colnames(datExpr) %in% choose_gene]
#view(datExpr)
#转置后存起来
datExpr2 <- as.data.frame(t(datExpr2))
write.csv(datExpr2,'yellow_moduel_final_expression_matrix_T_protein_coding_genes.csv',row.names = TRUE)
write.table(datExpr2,file="yellow_module_T_matrix.txt",sep="\t",quote = FALSE,col.names = T,row.names = T)
#不转置存起来
write.csv(datExpr2,'yellow_moduel_final_expression_matrix_protein_coding_genes.csv',row.names = TRUE)
####choose_gene存放的就是这三个与肿瘤显著相关的模块所包含的基因,将基因名保存到primal_gene_list.txt文件中

##计算TOM矩阵---这里只构建选中基因的TOM矩阵，代表它们之间的关联强度，以便后续作为基因调控网络的先验信息
adjacency = adjacency(datExpr2, power = power)  #这里使用datExpr2,不构建所有基因的TOM矩阵
TOM = TOMsimilarity(adjacency)
write.csv(TOM,'protein_coding_genes_TOM_Matrix.csv',row.names = TRUE)


#####################绘制模块之间相关性图#################################################
rm(list=ls())
load(file = "genes_modules.Rdata")
load(file = "pca_inspecting_after.Rdata")
load(file = "step4_design.Rdata")

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs   #每个模块的主成分基因
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

###########################计算模块和基因相关度，性状和基因相关度，画MS图###########################
### 计算模块与基因的相关性矩阵
corType="pearsoon"

if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(datExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}

# 计算性状与基因的相关性矩阵
##只有连续型性状才能进行计算，如果是离散变量，在构建样品表性时就转为0-1矩阵。

if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(datExpr, design, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "yellow"
pheno = "tumor"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design))
# 获取模块内的基因
moduleGenes = moduleColors == module
table(moduleGenes)
yellow_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(yellow_module)="genename"

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "grey",pch=20)
abline(h=0.5,v=0.8,col="red",lwd=1.5)

#按照指标选出hub基因集
module = "yellow"
column = match(module, modNames)
moduleGenes = moduleColors==module
lightyellow_module<-as.data.frame(dimnames(data.frame(datExpr))[[2]][moduleGenes])
names(lightyellow_module)="genename"
MM<-abs(geneModuleMembership[moduleGenes,column])
GS<-abs(geneTraitCor[moduleGenes, 1])

c<-as.data.frame(cbind(MM,GS))
rownames(c)=yellow_module$genename
head(c)

lightyellow_MMGS<-as.data.frame(cbind(MM,GS))
rownames(lightyellow_MMGS)=lightyellow_module$genename
hub_b<-abs(lightyellow_MMGS$MM)>0.8&abs(lightyellow_MMGS$GS)>0.2
table(hub_b)

lightyellow_hub_b<-subset(lightyellow_MMGS, abs(lightyellow_MMGS$MM)>0.8&abs(lightyellow_MMGS$GS)>0.2)
write.csv(lightyellow_hub_b, "hubgene_MMGS_lightyellow.csv")

yellow_hub<-read.csv("hubgene_MMGS_lightyellow.csv")
head(yellow_hub)
match<- yellow_module$genename %in% yellow_hub$X
c$group<-match
head(c)
library(ggplot2)
#pdf("MM vs. GS_yellow_TL.pdf",width = 7,height = 7)
ggplot(data=c, aes(x=MM, y=GS, color=group))+geom_point(size=1.5)+scale_colour_manual(values=c("grey60", "#DE6757"))+ theme_bw()+  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+  labs(x="Module Membership in yellow module", y="Gene significance for Tumor",title = "Module membership vs. Gene significance ")+theme(axis.title.x =element_text(size=14), axis.title.y=element_text(size=14),axis.text = element_text(size = 12),axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"),plot.title = element_text(hjust = 0.5,size = 16,face = "bold"),plot.margin = unit(rep(2,4),'lines')) +theme(legend.position = 'none')+geom_hline(aes(yintercept=0.2),colour="#5B9BD5",lwd=1,linetype=5)+geom_vline(aes(xintercept=0.8),colour="#5B9BD5",lwd=1,linetype=5)
#dev.off()