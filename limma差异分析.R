group <- c(rep("星形胶质细胞缺氧",3),rep("星形胶质细胞常氧",3),rep("HeLa缺氧",3),rep("HeLa常氧",3))







designtable<-model.matrix(~0+factor(group))
colnames(designtable)<-levels(factor(group))
rownames(designtable)<-colnames(exp_symbol)
#contrast matrix
contrast.matrix=makeContrasts(星形胶质细胞缺氧-星形胶质细胞常氧-HeLa缺氧-HeLa常氧,levels=designtable)
#limma DEG
fit<-lmFit(exp_symbol,designtable)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
options(digits=4)
DEG<-topTable(fit2,coef="星形胶质细胞缺氧 - 星形胶质细胞常氧 - HeLa缺氧 - HeLa常氧",n=Inf)#根据比较的组修改coef值
#DEG<-topTable(fit2,coef=3,n=Inf) 根据比较的组数修改coef值
DEG$group<-ifelse(DEG$P.Value>0.05,"no_change",
                  ifelse(DEG$logFC>1,"up",
                         ifelse(DEG$logFC< -1,"down","no_change")))
table(DEG$group)
DEG$gene<-rownames(DEG)
write.table(DEG,"DEG.txt",sep="\t",#写出文件
            quote=F,
            col.names=T,row.names = F)

