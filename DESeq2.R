##Check Opposum
library("DESeq2")
library("ggplot2")
library("dplyr")
library("ggrepel") #Avoid overlapping labels
library("plotly")
library("gplots")
library("VennDiagram")


Sys.setenv("plotly_username"="said3427")
Sys.setenv("plotly_api_key"="zVjG3VJ2tKDk1ed7FYhZ")


data<-read.table("merged_gene_name_counts.txt",header=T)
rownames(data)<-data[,1]
data<-data[,-1]
colnames(data)<-unlist(strsplit(colnames(data),"_gene_name.featureCounts.txt"))

data<-data[,-5]

#Reorder samples
data<-data[,c(9:11,1:8)]

group=c(rep("NL",3),rep("D12",4),"D18","D18R2CD2",rep("D18",2))
#group=c(rep("D12",4),rep("D18",5),rep("NL",3))

coldata<-data.frame(condition=group,type="paired-end")
rownames(coldata)<-colnames(data)

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design = ~ condition)




dds <- dds[ rowSums(counts(dds)) > 1, ]
dds$condition <- relevel(dds$condition, ref="NL")
pdf("boxplotPreNormalization.pdf")
boxplot(log2(counts(dds)+1),las=2,col=c(rep("green",3),rep("blue",4),rep("red",4)))
dev.off()
#dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)

conteos<-counts(dds,normalized=TRUE)


plotMA(dds,ylim=c(-7,7))


D18R_D12R.deseq2 <- results(dds, contrast=c("condition","D18","D12"),alpha=.05,independentFiltering=FALSE)
D18R_D12R.deseq2 <- D18R_D12R.deseq2[order(D18R_D12R.deseq2$pvalue),]

D18R_NLR.deseq2 <- results(dds, contrast=c("condition","D18","NL"),alpha=.05,independentFiltering=FALSE,pAdjustMethod="BH")
D18R_NLR.deseq2 <- D18R_NLR.deseq2[order(D18R_NLR.deseq2$pvalue),]

D12R_NLR.deseq2 <- results(dds, contrast=c("condition","D12","NL"),alpha=.05,independentFiltering=FALSE,pAdjustMethod="BH")
D12R_NLR.deseq2 <- D12R_NLR.deseq2[order(D18R_NLR.deseq2$pvalue),]


D18R_D12R <- data.frame(ID=rownames(D18R_D12R.deseq2),D18R_D12R.deseq2)
D18R_NLR <- data.frame(ID=rownames(D18R_NLR.deseq2),D18R_NLR.deseq2)
D12R_NLR <- data.frame(ID=rownames(D12R_NLR.deseq2),D12R_NLR.deseq2)

id.D18R_D12R<-subset(D18R_D12R,abs(log2FoldChange) > 1.5 & pvalue <.05)$ID
id.D18R_NLR<-subset(D18R_NLR,abs(log2FoldChange) > 1.5 & padj <.05)$ID
id.D12R_NLR<-subset(D12R_NLR,abs(log2FoldChange) > 1.5 & padj <.05)$ID



dim(subset(D18R_NLR,abs(log2FoldChange) > 4 & padj <.05))

write.csv(conteos,"conteos.csv",quote=FALSE)
save(conteos,id.D12R_NLR,id.D18R_NLR,id.D18R_D12R,id.union,file="heatmap.RData")

id.D12R_NLR<-id.D12R_NLR[-49]

#par(mar=c(8,4,1,2)+0.1) 

pdf("Heatmap_D18-D12.pdf",height = 12,width = 10);
heatmap.2(subset(conteos,rownames(conteos)%in%id.D18R_D12R),main ="D18-D12", ColSideColors = c(rep("green",3),rep("blue",4),rep("red",4)),scale="row", col=greenred(10),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none",Colv = NULL,dendrogram = "none",margins=c(6,10))
dev.off()

pdf("Heatmap_D18-NL.pdf",height = 22,width = 10);
heatmap.2(subset(conteos,rownames(conteos)%in%id.D18R_NLR),main="D18-NL",ColSideColors = c(rep("green",3),rep("blue",4),rep("red",4)), scale="row", col=greenred(10),symm=F,symkey=TRUE,symbreaks=FALSE,trace="none",Colv = NULL,dendrogram = "none",cexRow = .075, margins=c(6,10),offsetRow = -.3)
dev.off()

pdf("Heatmap_D12-NL.pdf",height = 22,width = 10);
heatmap.2(subset(conteos,rownames(conteos)%in%id.D12R_NLR),main="D12-NL",ColSideColors = c(rep("green",3),rep("blue",4),rep("red",4)), scale="row",col=greenred(10),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none",Colv = NULL, dendrogram = "none",cexRow = .075, margins=c(6,10),offsetRow = -.3)
dev.off()


heatmap.2(subset(conteos,rownames(conteos)%in%id.D12R_NLR),ColSideColors = c(rep("blue",3),rep("red",4),rep("green",4)), scale = "row",col=greenred(10),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none", Colv = NULL,dendrogram = 'none')

library(gplots)
pdf("GenesTopD18D12.pdf",height = 12,width = 10);
heatmap.2(subset(conteos,rownames(conteos)%in%id.union),ColSideColors = c(rep("green",3),rep("blue",4),rep("red",4)), scale = "row",col=greenred(100),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none",dendrogram = "row",Colv = "none",density.info = "none",margins = c(12,10)); 
dev.off()

pdf("HeatmapD18D12AllSamples.p
    df",height = 12,width = 10);
heatmap.2(subset(conteos,rownames(conteos)%in%id.D18R_D12R),ColSideColors = c(rep("green",3),rep("blue",4),rep("red",4)), scale = "row",col=greenred(100),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none",density.info = "none",margins = c(12,10))
dev.off()

pdf("HeatmapD18D12.pdf",height = 12,width = 10);
heatmap.2(subset(conteos[,-1:-3],rownames(conteos)%in%id.D18R_D12R),ColSideColors = c(rep("blue",4),rep("red",4)), scale = "row",col=greenred(100),symm=F,symkey=TRUE,symbreaks=TRUE,trace="none",density.info = "none",margins = c(12,10))
dev.off()




rld <- rlog(dds, blind=FALSE)
dists <- dist(t(assay(rld)))
plot(hclust(dists))


pca <- prcomp(t(assay(rld)))
plot(pca$x,pch=""); text(pca$x,labels = rownames(pca$x),cex=.5)

results<-D18R_D12R
results["Gene"]<-rownames(results)
results = mutate(results, sig=ifelse(results$pvalue<0.05 & log2FoldChange>1.5, "Sobreexpresados", ifelse(results$pvalue<0.05 & log2FoldChange< -1.5,"Subexpresados"," ")))
results = mutate(results, sig=ifelse((results$pvalue<0.05 & results$log2FoldChange>1.5), "Sobreexpresados", ifelse((results$pvalue<0.05 & results$log2FoldChange< -1.5),"Subexpresados","No Significativos")))
results$text=paste("Gen: ",results$Gene, "<br />log2(Tasa de Cambio): ", round(results$log2FoldChange,digits = 2),"<br />valor-p: ",formatC(results$pvalue,format="e",digits=2),sep="")

p <- plot_ly(results, x=~log2FoldChange,y=~-log10(pvalue),text=~text,color=~sig,type="scatter",mode="markers",colors=c("black","red","green")) %>% layout(xaxis=list(title="log2( Tasa de Cambio )"),yaxis=list(title="-log10(valor p)"),title="D18-D12")
p


p = ggplot(results, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col=sig)) +
  scale_color_manual(values=c("red", "black"))+ geom_text(aes(label=Gene),size=0)+
  ggtitle("D18-NL") 

p+theme_classic()


#p+geom_text_repel(data=filter(results, pvalue<0.05 & abs(log2FoldChange)>1.5), aes(label=Gene))

ggsave("VolcanoplotD18_NL.jpeg", device="jpeg") #In case you want to easily save to disk










p3d <- plot_ly(df.plot, x= ~PC1, y= ~PC2, z= ~PC3, color = ~class, colors = c("blue","green","magenta")) %>% add_markers()%>%  layout(scene = list(xaxis = list(title = 'PC1'),
                                                                                                                                                 yaxis = list(title = 'PC2'), zaxis = list(title = 'PC3')))


p2d <- plot_ly(df.plot, x= ~PC1, y= ~PC2, color = ~class, colors = c("blue","green","magenta")) %>% add_markers()%>%  layout(scene = list(xaxis = list(title = 'PC1'),
                                                                                                                                                   yaxis = list(title = 'PC2')))

venn.diagram(x=list(id.D18R_NLR,id.D12R_NLR,id.D18R_D12R),category.names = c("D18 - NL", "D12 - NL", "D18 - D12"), output=TRUE,filename = "Venn_poles.png",fill = c('red', 'blue', 'green'))

