library(corrplot)
library(ggplot2)
library(cowplot)
library(Hmisc)
library(pheatmap)
library(WGCNA)
library(RColorBrewer)
mycolor<-brewer.pal(9,"Set1")
options(bitmapType='cairo')
col <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(200)
ensg2gene <- read.table("/home/yuanjian/Project/ESCC_epigenome/Homeobox/ensg2gene")
gene_biotype <-
  read.table(
    "/home/yuanjian/Project/ESCC_epigenome/Homeobox//Homo_sapiens.GRCh38.95.genebiotype",
    sep = "\t",
    as.is = T,
    col.names = c("gene", "biotype")
  )
pcg <- gene_biotype[which(gene_biotype$biotype=="protein_coding"),]
pcg <- merge(pcg,ensg2gene,by.x="gene",by.y="V1")
#correlation----
# Cornea
Cornea.data <- read.delim("/home/yuanjian/Project/COVID/GSE77938.Corneas.final.TPM.coding.txt", header = T)
Cornea.data<- Cornea.data[,which(colnames(Cornea.data) %in% setdiff(colnames(Cornea.data),c("KR_19","KR_21","KR_24","KR_25","KR_26","KR_35")))]
dim(Cornea.data)
select1 <- t(Cornea.data)[,c("ACE2", "TMPRSS2", "SLC6A20", "LZTFL1", "FYCO1", "CCR9", "CXCR6", "XCR1")]
cor1 <- cor(select1)
res1 <- cor.mtest(select1, conf.level = .95)
corrplot(cor1, method = "circle", type = "upper", is.corr = T, diag = T, tl.col = "black", order = "original", tl.cex = 1,
	cl.cex = 1, cl.align.text = "l", p.mat = res1$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.cex = 1.5, col = col)
# Lung
Lung.data <- read.delim("GTEx.Lung.final.TPM.txt", header = T)
dim(Lung.data)
select2 <- t(Lung.data)[,c("ACE2", "TMPRSS2", "SLC6A20", "LZTFL1", "FYCO1", "CCR9", "CXCR6", "XCR1")]
cor2 <- cor(select2)
res2 <- cor.mtest(select2, conf.level = .95)
corrplot(cor2, method = "circle", type = "upper", is.corr = T, diag = T, tl.col = "black", order = "original", tl.cex = 1,
         cl.cex = 1, cl.align.text = "l", p.mat = res2$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.cex = 1.5, col = col)
# RPE
RPE.data <- read.delim("RPE.final.nRPKM.txt", header = T)
dim(RPE.data)
select3 <- t(RPE.data)[,c("ACE2", "TMPRSS2", "SLC6A20", "LZTFL1", "FYCO1", "CCR9", "CXCR6", "XCR1")]
cor3 <- cor(select3)
res3 <- cor.mtest(select3, conf.level = .95)
corrplot(cor3, method = "circle", type = "upper", is.corr = T, diag = T, tl.col = "black", order = "original", tl.cex = 1,
         cl.cex = 1, cl.align.text = "l", p.mat = res3$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", pch.cex = 1.5, col = col)
#heatmap----
GTEx <- read.delim("GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", header = T, skip = 2,row.names = 1)
head(GTEx[,1:5])
tissue <- read.delim("GTEx.tissues.txt", header = T, skip = 1)
gene <- c("ACE2", "TMPRSS2", "SLC6A20", "LZTFL1", "FYCO1", "CCR9", "CXCR6", "XCR1")
data <- cbind(GTEx[match(gene, GTEx$Description),-1], apply(Cornea.data[match(gene, rownames(Cornea.data)),], 1, median),
              apply(RPE.data[match(gene, rownames(RPE.data)),], 1, median))
colnames(data) <- c(as.character(tissue[,1]), "Cornea", "RPE")
rownames(data) <- gene
data <- t(t(data)[order(t(data)[,"ACE2"], decreasing = F),])
heatmap.data <- log2(data+1)
col <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"))(200)
par(mar = c(6, 7, 5, 2), lwd = 2)
pheatmap(heatmap.data, cluster_rows = F, cluster_cols = F, show_rownames = T, show_colnames = T,
         border_color = NA, scale = "none", fontsize = 13, col = col, angle_col = "45", width = 16, height = 5)

#kmeans clustering----
kmeans <- kmeans(data,opt$cluster,iter.max=1000,nstar=10)
kmeans <- read.table("/home/yuanjian/Project/COVID/kmeans_cluster.txt",header = T)
datExpr0_combine <- read.table("/home/yuanjian/Project/COVID/GSE77938.Corneas.final.TPM.coding.txt")
datExpr0_combine <- datExpr0_combine[which(rownames(datExpr0_combine) %in% kmeans$id),]
datExpr0_combine <- datExpr0_combine[which(rownames(datExpr0_combine) %in% pcg$V2),]
datExpr0_combine<- datExpr0_combine[,which(colnames(datExpr0_combine) %in% setdiff(colnames(datExpr0_combine),c("KR_19","KR_21","KR_24","KR_25","KR_26","KR_35")))]
percentile <- function(x) {
  length(x[x >= 1]) / 19
}
datExpr0_combine$percentile <- apply(datExpr0_combine, 1, percentile)
datExpr0_combine <-
  datExpr0_combine[which(datExpr0_combine$percentile >= 0.05), 1:19]
datExpr1_combine <- as.matrix(datExpr0_combine)
datExpr1_combine <- normalize.quantiles(datExpr1_combine,copy = TRUE)
colnames(datExpr1_combine) <- colnames(datExpr0_combine)
rownames(datExpr1_combine) <- rownames(datExpr0_combine)
datExpr2_combine <- data.frame(t(datExpr0_combine))
colnames(datExpr2_combine) <- rownames(datExpr1_combine)
dim(datExpr2_combine)
#[1] 19 18274
datExpr3_combine <- data.frame(t(datExpr0_combine))
colnames(datExpr3_combine) <- rownames(datExpr0_combine)
dim(datExpr3_combine)
#[1] 19 18274
all_gene <- ensg2gene[which(ensg2gene$V2 %in% colnames(datExpr2_combine)),]
write.table(all_gene$V1,file = "all.txt",quote = F,col.names = F,row.names = F)
kmeans <- kmeans[which(kmeans$id %in% colnames(datExpr2_combine)),]
kmeans <- kmeans[order(match(kmeans$id,colnames(datExpr2_combine))),]
dim(kmeans)
#[1] 18274    21
moduleLabels_combine_signed = kmeans$cluster
moduleColors_combine_signed = labels2colors(moduleLabels_combine_signed)
length(table(moduleColors_combine_signed))
gene_color <- as.data.frame(cbind(colnames(datExpr2_combine),moduleLabels_combine_signed,moduleColors_combine_signed))
table(gene_color$moduleLabels_combine_signed)
gene_color_ace2 <- gene_color[which(gene_color$moduleLabels_combine_signed==7),]
cluster_gene <- ensg2gene[which(ensg2gene$V2 %in% gene_color_ace2$V1),]
write.table(cluster_gene$V1,file = "cluster_gene.txt",quote = F,col.names = F,row.names = F)
gene_color[which(gene_color$V1 %in% c("ACE2","TMPRSS2","SLC6A20","LZTFL1","FYCO1","CCR9","CXCR6","XCR1")),]
MEs0 = moduleEigengenes(datExpr2_combine,
                        moduleColors_combine_signed,
                        nPC = 1,
                        excludeGrey = T)$eigengenes
cor1 <- cor(MEs0)
res1 <- cor.mtest(MEs0, conf.level = .95)
corrplot(cor1, method = "square", is.corr = T, diag = T, tl.col = "black", order = "hclust", tl.cex = 1,
         cl.cex = 1, cl.align.text = "l", insig = "blank", pch.cex = 1.5, col = col)
MEs = orderMEs(MEs0)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr2_combine, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
geneModuleMembership[which(rownames(geneModuleMembership) %in% c("ACE2","TMPRSS2","SLC6A20","LZTFL1","FYCO1","CCR9","CXCR6","XCR1")),]
geneModuleMembership[which(rownames(geneModuleMembership) %in% c("MRPS5","MRPS27","MRPS25","MRPS2")),]

geneModuleMembership_ACE2 <- geneModuleMembership[which(rownames(geneModuleMembership) %in% gene_color_ace2$V1),]
geneModuleMembership <- geneModuleMembership[order(geneModuleMembership$MMblack,decreasing = T),]
geneModuleMembership <- geneModuleMembership[order(geneModuleMembership$MMgreen,decreasing = T),]
geneModuleMembership <- geneModuleMembership[order(geneModuleMembership$MMbrown,decreasing = T),]

a <- geneModuleMembership[which(geneModuleMembership$MMblack >= 0.5),]
a$max <- apply(a,1,max)
a <- a[which(a$max==a$MMblack),]
b <- geneModuleMembership_ACE2[which(rownames(geneModuleMembership_ACE2) %in% rownames(a)),]
b[which(rownames(b) %in% c("ACE2","TMPRSS2","SLC6A20","LZTFL1","FYCO1","CCR9","CXCR6","XCR1")),]
cluster_gene1 <- cluster_gene[which(cluster_gene$V2 %in% rownames(b)),]
write.table(cluster_gene1$V1,file = "cluster_gene1.txt",col.names = F,row.names = F,quote = F,sep = "\t")

length(intersect(rownames(a),mito_gene$V2))
length(intersect(rownames(b),mito_gene$V2))
length(intersect(cluster_gene$V2,mito_gene$V2))

length(intersect(rownames(a),inter_gene$V2))
length(intersect(rownames(b),inter_gene$V2))
length(intersect(cluster_gene$V2,inter_gene$V2))


interactome <- data.frame(gene=kmeans$id,cluster=1)
interactome1 <- interactome[which(interactome$gene %in% inter$V1),]
interactome1$cluster <- 2
interactome2 <- interactome[which(interactome$gene %in% setdiff(interactome$gene,interactome1$gene)),]
interactome3 <- rbind(interactome1,interactome2)
moduleLabels_combine_signed = interactome3$cluster
moduleColors_combine_signed = labels2colors(moduleLabels_combine_signed)
length(table(moduleColors_combine_signed))
MEs1 = moduleEigengenes(datExpr2_combine,
                        moduleColors_combine_signed,
                        nPC = 1,
                        excludeGrey = T)$eigengenes

cor.test(MEs0$MEblack,MEs1$MEblue)



interactome <- data.frame(gene=kmeans$id,cluster=1)
candidate <- intersect(cluster_gene$V2,mito_gene$V2)
candidate <- intersect(candidate,inter_gene$V2)
interactome1 <- interactome[which(interactome$gene %in% candidate),]
interactome1$cluster <- 2
interactome2 <- interactome[which(interactome$gene %in% setdiff(interactome$gene,interactome1$gene)),]
interactome3 <- rbind(interactome1,interactome2)
moduleLabels_combine_signed = interactome3$cluster
moduleColors_combine_signed = labels2colors(moduleLabels_combine_signed)
length(table(moduleColors_combine_signed))
MEs1 = moduleEigengenes(datExpr2_combine,
                        moduleColors_combine_signed,
                        nPC = 1,
                        excludeGrey = T)$eigengenes
MEs = orderMEs(MEs1)
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr2_combine, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")
geneModuleMembership[which(rownames(geneModuleMembership) %in% candidate),]
#GO correlation----
inter <- read.table("/home/yuanjian/Project/COVID/inter.txt")
inter_gene <- ensg2gene[which(ensg2gene$V2 %in% inter$V1),]
write.table(inter_gene$V1,file = "inter_gene.txt",quote = F,col.names = F,row.names = F)

inter_GO <- read.table("/home/yuanjian/Project/COVID/inter_gene_GO.txt",sep = "\t",header = T)
cluster_GO <- read.table("/home/yuanjian/Project/COVID/cluster_gene_GO.txt",sep = "\t",header = T)
merge_GO <- merge(inter_GO,cluster_GO,by="Term")
merge_GO <- merge_GO[order(merge_GO$Fold.Enrichment.x),]
par(mar = c(6, 7, 5, 2), lwd = 3)
ggplot(merge_GO, aes(x=Fold.Enrichment.x, y=Fold.Enrichment.y)) + 
  geom_point(color='#2980B9', size = 4) + 
  geom_smooth(method=lm, color='#2C3E50')


#venn----
mito <- read.csv("/home/yuanjian/Project/COVID/results.csv",header = F)
mito_gene <- ensg2gene[which(ensg2gene$V1 %in% mito$V1),]

venn <- list(
  a=cluster_gene$V1,
  b=mito_gene$V1,
  c=inter_gene$V1)
names(venn)<-c("cluster","mito","inter")
venn.diagram(venn,filename="venn.pdf",height=1000,width=1000,fill=mycolor[1:3],cex=0.8,fontface=2,lty=0)


length(intersect(cluster_gene$V1,mito_gene$V1))
random_cluster <- sample(rownames(datExpr0_combine),1435)
random_cluster <- ensg2gene[which(ensg2gene$V2 %in% random_cluster),]
length(intersect(random_cluster$V1,mito_gene$V1))
fisher.test(matrix(c(221,1435,83,1435),nrow = 2))


length(intersect(cluster_gene$V1,inter_gene$V1))
random_cluster <- sample(rownames(datExpr0_combine),1435)
random_cluster <- ensg2gene[which(ensg2gene$V2 %in% random_cluster),]
length(intersect(random_cluster$V1,inter_gene$V1))
fisher.test(matrix(c(41,1435,20,1435),nrow = 2))

length(intersect(mito_gene$V1,inter_gene$V1))
random_cluster <- sample(rownames(datExpr0_combine),1157)
random_cluster <- ensg2gene[which(ensg2gene$V2 %in% random_cluster),]
length(intersect(random_cluster$V1,inter_gene$V1))
fisher.test(matrix(c(56,1157,16,1157),nrow = 2))

x<-ensg2gene[which(ensg2gene$V2 %in% intersect(intersect(mito_gene$V2,inter_gene$V2),cluster_gene$V2)),]

#drug----
drug <- read.csv("drugbank_mitochondrial.csv")
drug1 <- read.csv("drug1.csv")
drug2 <- read.csv("drug2.csv")

candidate <- intersect(intersect(cluster_gene$V2,mito_gene$V2),inter_gene$V2)
write.table(candidate,file = "candidate.txt",col.names = F,row.names = F,quote = F,sep = "\t")
drug[which(drug$Gene.Name %in% candidate),]
drug1[which(drug1$prey_gene_name %in% candidate),]
drug2[which(drug2$prey_gene_name %in% candidate),]


candidate1 <- intersect(cluster_gene$V2,mito_gene$V2)
a1 <- drug[which(drug$Gene.Name %in% candidate1),]
a2 <- drug1[which(drug1$prey_gene_name %in% candidate1),]
a3 <- drug2[which(drug2$prey_gene_name %in% candidate1),]

target <- read.csv("target.csv")
target_c <- target[which(target$Saint_BFDR <= 0.05 & target$AvgSpec>=2 & target$MIST>=0.7),]

#DEG----
a <- read.table("Project/COVID/Keratoconus_vs_Cornea.DEG.TPM.txt")
a1 <- a[which(a$Group =="Corneas"),]
a2 <- a[which(a$Group =="Keratoconus"),]
t.test(a1$ACE2,a2$ACE2)
ks.test(a1$TMPRSS2,a2$TMPRSS2)
ks.test(a1$SLC6A20,a2$SLC6A20)
t.test(a1$CCR9,a2$CCR9)
wilcox.test(a1$FAM162A,a2$FAM162A)
wilcox.test(a1$DNAJC19,a2$DNAJC19)
t.test(a1$SLC30A9,a2$SLC30A9)
ks.test(a1$NDUFB9,a2$NDUFB9)
t.test(a1$NDUFAF1,a2$NDUFAF1)
t.test(a1$ECSIT,a2$ECSIT)