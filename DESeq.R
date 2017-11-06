## Clean the working space ---------------------------------------------
rm (list = ls ())

## Install packages ----------------------------------------------------
source("http://www.bioconductor.org/biocLite.R")
#biocLite("ALL")
biocLite("DESeq")
biocLite("GOstats")
#biocLite("DESeq2")

## Load packages -------------------------------------------------------
library(DESeq)
library(GOstats)
library(GO.db)
if (!require("gplots")) { install.packages("gplots", dependencies = TRUE) }
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) }
if (!require("reshape2")) { install.packages("reshape2", dependencies = TRUE) }
#library("DESeq2")

## Import Data ---------------------------------------------------------
setwd("~/Data/GeneFams/ESPs/HM-Expression/")
countsTable <- read.delim("Counts.txt",header=TRUE)
rownames(countsTable) <- countsTable$Gene
countsTable <- countsTable[,-1]
genotype <- factor(c(rep("HM034mel",3), rep("HM056mel",3), rep("HM101mel",3), rep("HM340mel",3),rep("HM034med",3), rep("HM056med",3), rep("HM101med",3), rep("HM340med",3)))
cds <- newCountDataSet( countsTable, genotype)

## Normalize the data -------------------------------------------------
cds <- estimateSizeFactors(cds)
sizeFactors(cds)
cds <- estimateDispersions(cds)
norm <- counts(cds, normalized=TRUE)

##
countsTable[14,]/13255161
mysizefactors = 
sizeFactors(cds) = c(my_Values)
#x <- sizeFactors(cds)
#sizeFactors(cds2) <- x
#cds <- estimateDispersions(cds2)
#norm <- counts(cds2, normalized=TRUE)



#####Test for DE genes --------------------------------------------------------
#Groups: nodHM nod12 nod85 nod26
#test1=3,5=Gm12.1   Test2=3,5,6=Gm8.5   Test3=2,3,5=Gm26.1
test1 = nbinomTest(cds, "HM034mel", "HM034med")
test2 = nbinomTest(cds, "HM056mel", "HM056med")
test3 = nbinomTest(cds, "HM101mel", "HM101med")
test4 = nbinomTest(cds, "HM340mel", "HM340med")
format(test1, scientific=FALSE)
format(test2, scientific=FALSE)
format(test3, scientific=FALSE)
format(test4, scientific=FALSE)
#subset(test5, padj<0.1 & !(foldChange >= 0.5 & foldChange <= 2))
subset(test4, id=="Medtr4g094732")


# Heatmap of all DE genes with padj<0.1 ---------------------------------------
#AllDE = scan("~/Data/PDPs/DESeq/DElist.txt", what="", sep="\n")
AllDE = rownames(norm)[c(1:11,13)]
DEtable = norm[row.names(norm) %in% AllDE, ]
ordercols <- colnames(DEtable)[c(1,2,3,13,14,15,4,5,6,16,17,18,7,8,9,19,20,21,10,11,12,22,23,24)]
DEtable = DEtable[,match(ordercols, colnames(DEtable))]
ScaleDE = scale(t(DEtable))
#Yellow indicates high  expression while red indicates low expression
heatmap.2(t(ScaleDE), dendrogram = 'none', Rowv = NULL, Colv=NULL, margin=c(8,8), key=FALSE,  trace='none')


## Plot correlations of normalized expression levels ------------------
cor<-cor(t(ScaleDE), use="pairwise.complete.obs", method="pearson")
#pdf("heatmap.pdf", pointsize=10, font="Helvetica")
heatmap.2(cor, trace="none", col="redgreen", density.info="none", Rowv=FALSE, Colv=FALSE, dendrogram="none", margin=c(6,6))
heatmap.2(cor, trace="none", col="redgreen", density.info="none", margin=c(6,6))
#dev.off()
hist(cor, breaks=50)
dist<-as.dist((1-cor)/2)
#pdf(file="tree.pdf")
plot(hclust(dist))
#dev.off()


# Interaction plots --------------------------------------------------
df <- melt(DEtable, varnames=c("Gene","Sample"))
genotype <- factor(c(rep("HM034",72), rep("HM056",72), rep("HM101",72), rep("HM340",72)))
rhizobia <- factor(c(rep(c(rep("S.meliloti",36),rep("S.medicae",36)),4)))
df<-cbind(genotype, rhizobia, df)

#Plot
subsetdf <- subset(df, df$Gene == "PDP3")
interaction.plot(subsetdf$genotype, subsetdf$rhizobia, subsetdf$value, col = "purple", 
                 trace.label = "Strain", fixed = TRUE, xlab = "Genotype", ylab = "Counts",
                 type = "b",pch=c(21,19))

se=function(x) sqrt(var(x)/length(x)) 
#Calculate the means and standard error for each type x treatment 
subsetdf <- subset(df, df$Gene == "PDP5")
label = "NPD4"
means=by(subsetdf$value,list(subsetdf$rhizobia,subsetdf$genotype),mean) 
ses=by(subsetdf$value,list(subsetdf$rhizobia,subsetdf$genotype),se) 
sd=by(subsetdf$value,list(subsetdf$rhizobia,subsetdf$genotype),sd) 
#plot means as interaction plot; type="b" means plot both symbols and lines; pch=c(21,19) are the two symbol types; ylim is the y axis minimum & maximum values; las=1 makes the y-axis numbers horizontal 
interaction.plot(subsetdf$genotype,subsetdf$rhizobia,subsetdf$value, col = "purple",legend=TRUE, 
                 trace.label = "Rhizobial strain", fixed = TRUE, xlab = "Medicago Accession", ylab = "Counts",
                 type="b",pch=c(21,19), lwd = 1.5, cex = 1.2, 
                 ylim = c(min(subsetdf$value)+50, max(subsetdf$value)-50)) 
#now add the standard error lines (means plus/minus standard error) 
title(label, line = -9, adj=0.9, col.main = "seagreen4", cex.main = 2.2)
lines(c(1,1),c(means[1]-ses[1],means[1]+ses[1]), col = "purple", lwd = 1.5) 
lines(c(1,1),c(means[2]-ses[2],means[2]+ses[2]), col = "purple", lwd = 1.5) 
lines(c(2,2),c(means[3]-ses[3],means[3]+ses[3]), col = "purple", lwd = 1.5) 
lines(c(2,2),c(means[4]-ses[4],means[4]+ses[4]), col = "purple", lwd = 1.5)
lines(c(3,3),c(means[5]-ses[5],means[5]+ses[5]), col = "purple", lwd = 1.5) 
lines(c(3,3),c(means[6]-ses[6],means[6]+ses[6]), col = "purple", lwd = 1.5) 
lines(c(4,4),c(means[7]-ses[7],means[7]+ses[7]), col = "purple", lwd = 1.5) 
lines(c(4,4),c(means[8]-ses[8],means[8]+ses[8]), col = "purple", lwd = 1.5)







as.matrix(df)
# Heatmap of HM-Gm85(356) DE genes with padj<0.1 -----------------------------------
#DEgenes <- read.delim("~/Data/PDPs/DESeq/DEpdp356genes.txt", header=FALSE)
DEgenes = Sub.Ad2
DEtable = norm[row.names(norm) %in% DEgenes, ]
DEtable = DEtable[, c(1:9)]
ScaleDE = scale(t(DEtable))
heatmap(t(ScaleDE), Colv=FALSE)
#Save heatmap
png(filename='~/Data/PDPs/DESeq/DEpdp356-1.png', width=1800, height=1000)
heatmap.2(t(ScaleDE),dendrogram='row', Rowv=TRUE, Colv=FALSE, margin=c(12,80),
          cexRow = 2.8, cexCol = 2.5, key=FALSE, trace='none',density.info="none")
graphics.off()
#Get Heatmap order
Annotations = read.delim("~/Data/Mt4.0/Data/genome/MtDescriptionList.txt", header=FALSE)
hm = heatmap.2(t(ScaleDE),dendrogram='row', Rowv=TRUE, Colv=(FALSE))
hclust = as.hclust(hm$rowDendrogram)
geneorder = cutree(hclust, h=10)[hclust$order]
index = match(labels(geneorder), Annotations$V1)
labels = Annotations$V2[rev(index)]
#Labels
x <- 1:length(labels(geneorder))
y <- length(labels(geneorder)):1
png(filename='~/Data/PDPs/DESeq/DEpdp356-2.png', width=1800, height=1035)
 plot( x = x, y = y, pch = '', xlab = 'X',  ylab = 'Y' )
 text(12,y = -0.5+y, labels= labels, pos=4, cex=2.8)
graphics.off()

# Heatmap of HM-Gm12(pdp235) DE genes with padj<0.1 ---------------------------
#DEgenes <- read.delim("~/Data/PDPs/DESeq/DEpdp235.txt", header=FALSE)
DEgenes = Sub.Ad3
DEtable = norm[row.names(norm) %in% DEgenes$V1, ]
DEtable = DEtable[, c(1:6,10:12)]
ScaleDE = scale(t(DEtable))
#Save heatmap
png(filename='~/Data/PDPs/DESeq/DEpdp235-1.png', width=1800, height=4000)
heatmap.2(t(ScaleDE),dendrogram='row', Rowv=TRUE, 
          Colv=FALSE, margin=c(12,80),cexRow = 2.8, 
          cexCol = 2.5, key=FALSE, trace='none',density.info="none")
graphics.off()
#Get Heatmap order
hm = heatmap.2(t(ScaleDE),dendrogram='row', Rowv=TRUE, Colv=FALSE)
hclust = as.hclust(hm$rowDendrogram)
geneorder = cutree(hclust, h=10)[hclust$order]
index = match(labels(geneorder), Annotations$V1)
labels = Annotations$V2[rev(index)]
#Labels
y <- length(labels(geneorder)):1
png(filename='~/Data/PDPs/DESeq/DE235-3.png', width=1000, height=4100)
plot( x = c(1,length(labels(geneorder))), y = c(1,length(labels(geneorder))+10), pch = '', xlab = 'X',  ylab = 'Y' )
text(10,y = -3+y*1.04, labels=labels, pos=4, cex=2.8)
graphics.off()


##### Gene Ontology Enrichment  ---------------------------------------------
# GO enrichment of all DE genes with padj<0.1 -------------------------------
#AllDE = scan("~/Data/PDPs/DESeq/DEpdp235.txt", what="", sep="\n")
#test1=3,5=Gm12   Test2=3,5,6=Gm85   Test3=2,3,5=Gm26   Test4=12vs85   Test5=12vs26
setwd("~/Data/PDPs/GO")
source('~/Data/PDPs/GO/perform_enrichment.R')

AllDE = sort(Reduce(union, list(Sub.Pv1,Sub.Pv2,Sub.Pv3,Sub.Pv4,Sub.Pv5)))
write(AllDE, file = "AllDE.txt", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv05.1, file = "Gm12DE.txt", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv05.2, file = "Gm85DE.txt", ncolumns = 1, append = FALSE, sep = "\n")
write(Sub.Pv05.3, file = "Gm26DE.txt", ncolumns = 1, append = FALSE, sep = "\n")

TestGenes = scan("Gm26DE.txt", what="", sep="\n")

#Matrix backbone
AllMt = scan("~/Data/Mt4.0/Data/genome/UniqueList.txt", what="", sep="\n")
GeneMAT = matrix(data = 0, nrow = length(AllMt), ncol = 1)
row.names(GeneMAT) = c(AllMt)
colnames(GeneMAT) = c("DE")

# Put DE gene list into a matrix format. 1 is differentially expressed
for (i in 1:length(TestGenes)){
  temp <- which(TestGenes[i] == row.names(GeneMAT))
  GeneMAT[temp,1] <- 1
}


#######
##GOF is medicago go term propagation file
GOF <- read.delim("~/Data/PDPs/GO/GOMtMapped.txt", sep = '\t', stringsAsFactors = FALSE, header = FALSE)

colnames(GOF) <- c("Gene", "GO")
GOF <- cbind(GOF, value = 1)

xx <- as.list(GOBPANCESTOR)
xx <- xx[!is.na(xx)]

BIOP <- c(0)
for( i in 1:length(xx)){
  BIOP <- append(BIOP, names(xx[i]))
}

#restrict GO terms that are only in a biological process
final_mat = acast(GOF, Gene ~ GO, fill=0, value.var="value")


#filter for genes only from the ones tested
test = final_mat[,colnames(final_mat) %in% BIOP,drop=F]

#command to run the enrichments
Results <- perform_enrichment(GeneMAT, test, minSize = 3, p.value.cutoff=0.05)

#add function column on the end
Results$Function=Term(as.vector(Results$id))
write.table(Results, "Gm26enrichment.txt", quote = FALSE, sep='\t')


# Look at the number and plot DE genes --------------------------------------
test=test1
nrow(subset(test, pval<0.01 & !(foldChange >= 0.5 & foldChange <= 2)))
plot(test$log2FoldChange,-log10(test$pval))
plot(test$log2FoldChange~log10(test$baseMean), xlab="log10 Mean Gene Expression", ylab = "log2 Fold Change", main="HM340 vs. pdp3,pdp5")
#plot(test$log2FoldChange~log10(test$baseMean), xlab="log10 Mean Gene Expression", ylab = "log2 Fold Change", main="HM340 vs. pdp3,pdp5,pdp6")
points(test$log2FoldChange[test$pval<0.01]~log10(test$baseMean[test$pval<0.01]), col="red", pch=16)
points(test$log2FoldChange[test$padj<0.01]~log10(test$baseMean[test$padj<0.01]), col="blue", pch=16)     

# Choose lists of DE genes --------------------------------------------------
#test1=3,5=Gm12.1   Test2=3,5,6=Gm8.5   Test3=2,3,5=Gm26.1
#Padj < 0.1
Sub.Ad1 = subset(test1, padj<0.1 & !(foldChange >= 0.5 & foldChange <= 2))$id
#Pval < 0.01
Sub.Pv1 = subset(test1, pval<0.01 & !(foldChange >= 0.5 & foldChange <= 2))$id
#Pval < 0.05
Sub.Pv05.1 = subset(test1, pval<0.05 & !(foldChange >= 0.5 & foldChange <= 2))$id
