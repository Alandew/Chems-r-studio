

rm(list=ls())
setwd("C:\\Users\\湘灵\\Desktop\\dates2")
mat="GSE53757_series_matrix.txt"#基因表达矩阵
gpl="GPL570-55999.txt"#平台文件
a=1#平台文件的ID列号
b=11#平台文件的SYMBOL列号
gene_name="DDR1"#目标基因SYMBOL
out_name="data_out.txt"#输出文件名（记得加拓展名）
###################################################

#install.packages("pacman") 若未安装pacman包则执行此命令
library(pacman)
#library(GEOquery)加载数据包
p_load(limma,affy,tidyverse)
exp=read.table(mat,header=T,
               sep="\t",dec=".",
               comment.char="!",na.strings =c("NA"),fill=T )#读取基因表达矩阵
GPL_file=read.table(gpl,header=T,
                    quote="",sep="\t",dec=".",
                    comment.char="#",na.strings =c("NA"),fill=T )#读取平台文件（GPL）
gpl_file=GPL_file[,c(a,b)]#将平台文件的ID列和SYMBOL列取出
e <- apply(gpl_file,1,
           function(x){
             paste(x[1],
                   str_split(x[2],'///',simplify=T),
                   sep = "...")
           })
x = tibble(unlist(e))
colnames(x) <- "f" 
gpl_file <- separate(x,f,c("ID","symbol"),sep = "\\...")#处理一个探针对应多个基因
exp<-as.data.frame(exp)#将表达矩阵转换为数据框
colnames(exp)[1]="ID"#将表达矩阵的第一列列名改为ID（将表达矩阵和GPL的列名统一）
exp_symbol<-merge(exp,gpl_file,by="ID")#以ID为参照值，对表达矩阵和GPL进行合并
exp_symbol[exp_symbol==""]<-NA#将空白负值NA
exp_symbol<-na.omit(exp_symbol)#删除GENE_SYBOL缺失的数据
exp_symbol[,grep("symbol", colnames(exp_symbol))]=
  trimws(exp_symbol[,grep("symbol", colnames(exp_symbol))])#去除数据头尾空格
if("DDR1"%in% exp_symbol[,grep("symbol", colnames(exp_symbol))]== F)
{ errorCondition("未找到指定基因",class = NULL,call = NULL)}#检测是否含有目标基因
table(duplicated(exp_symbol[,ncol(exp_symbol)]))#检测重复基因数
d=data.frame(duplicated(exp_symbol[,ncol(exp_symbol)]))#将布尔值写入数据框
exp_symbol=avereps(exp_symbol[,-c(1,ncol(exp_symbol))],ID=exp_symbol$symbol)#对重复值取平均后混合
table(duplicated(rownames(exp_symbol)))#再次检测重复基因数
write.table(exp_symbol,out_name,sep="\t",#写出文件
            quote=F,
            col.names=T)
save(exp,exp_symbol,gpl_file,GPL_file,file="save_file.Rdata")
save(exp_symbol,file = 'exp_symbol.txt')
library(limma)
exprSet=normalizeBetweenArrays(exp_symbol)
qx <- as.numeric(quantile(exprSet, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)


## 开始判断
if (LogC) { 
  exprSet [which(exprSet  <= 0)] <- NaN
  ## 取log2
  exprSet_clean <- log2(exprSet)
  print("log2 transform finished")
}else{
  print("log2 transform not needed")
}


library("impute")


library(samr)

# 假设你的表达矩阵和分组信息已经准备好
expression_matrix <- as.matrix(exprSet_clean)
group_labels <- c(rep(c('unnormal','normal'),72))


library(siggenes)
samResult <- sam(expression_matrix, cl = group_labels)
significantGenes <- which(samResult@p.value < 0.05)
findDelta(fdr=0.05,samResult)
pzhi <- as.data.frame(significantGenes)
gene_name <- as.data.frame(rownames(pzhi))
colnames(gene_name) <- 'SAM'
write.csv(gene_name,file = 'SAM.txt')
pvalue <- as.data.frame(samResult@p.value)
fc <- as.data.frame(log2(samResult@fold))
nuew_file <- cbind(pvalue,fc)
colnames(nuew_file) <- c("p_value","logFC")
genname <- rownames(nuew_file[1:24428,])
nuew_file$gene <- genname
nuew_file$logFC <- logFC_matrix
t <- nuew_file[(nuew_file$p_value<0.05&(abs(nuew_file$logFC)>1)),]
up <- t[t$logFC>1,]
write.csv(t[1:1922,3],"sam1111.csv")


SampleInfo<-read.csv("sl.csv",header=F)#
SampleInfo$V1 <- NULL
colnames(SampleInfo) <- c(1:144)

logFC_matrix <- matrix(NA, nrow = nrow(exprSet_clean), ncol = 1)

# 计算logFC值
for(g in 1:(nrow(exprSet_clean))){
  k <- mean(exprSet_clean[g,which(SampleInfo[1,]=="tissue: clear cell renal cell carcinoma")])
  j <- mean(exprSet_clean[g,which(SampleInfo[1,]=="tissue: normal kidney")])
  logFC_matrix[g,1] <- (k-j)
  
  
  
}
write.csv(logFC_matrix,file = "logFC.csv")
# 现在，logFC_matrix 包含了所有基因针对每一对组的logFC值
p_value <- samResult@p.value
p_value <- as.data.frame(p_value)
p_value$logFC <- logFC_matrix[,1]
L <- p_value

colnames(L) <- c('pvalue','logFC')
data <- L
data$pvalue <- data$pvalue + 1e-300
data$gene <- rownames(data[1:3619,])
data <- data[(data$pvalue< 0.05 & (abs(data$logFC)>1)),]
data_up <- data[(data$logFC>1),]
data_down <- data[(data$logFC< -1),]

L$group<-ifelse(L$p_value>0.05,"no_change",
                  ifelse(L$logFC>1,"up",
                         ifelse(L$logFC< -1,"down","no_change")))
mm <- which(L$group=='up' | L$group=="down")
ppp <- which(L$group=='up')
m <- rownames(L[mm,])
m <- as.data.frame(m)
write.csv(m,"SAM包.csv")
# 假设你已经有了一个数据框，包含每个基因的logFC和-pvalue
# 这里我们构造一个示例数据框
library(ggplot2)
library(ggrepel)


# 合并上调和下调基因的数据  

# 分组信息
data <- t 


set.seed(123)
core_gene <- 0
core_gene
# 出图
ggplot(t,aes(logFC, -log10(p_value)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point(aes(size = -log10(p_value), 
                 color = -log10(p_value)))+
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"))+
  scale_size_continuous(range = c(0,1))+
  theme_bw(base_size = 12)+
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        legend.justification = c(0,1))+
  # 设置图例
  guides(col = 
           guide_colorbar(title = "-Log10_q-value",
                          ticks.colour = NA,
                          reverse = T,
                          title.vjust = 0.8,
                          barheight = 8,
                          barwidth = 1),
         size = "none") +
  # 添加标签：
  
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")# 绘制火山图  

