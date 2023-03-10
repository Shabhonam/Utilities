setwd("~/Desktop/Testdata/")
sampleInfo <- read.csv("experiment.csv", header=TRUE, row.names=1)
head(sampleInfo)
## add X at the beginning of rows beginning with a number (makes it consistent to column names of of the count matrix!)
if ( any(grepl("^[0-9]", sampleInfo$name)) ) {
  sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name = paste("X", sampleInfo[grepl("^[0-9]", sampleInfo$name),]$name, sep="")  
}
sampleInfo = DataFrame(as.data.frame(unclass(sampleInfo)))
##sampleInfo = sampleInfo[order(sampleInfo$name, decreasing=F),]  # order by sample name
as.character(sampleInfo$name)
head(sampleInfo)
## count matrix (e.g. from DESeq or featureCounts)

countdata<-read.csv("counts.csv",header=TRUE,row.names=1)
head(countdata)
countdata = DataFrame(countdata)
countdata = countdata[,as.character(sampleInfo[,1])]
head(countdata)
colnames(countdata)

dds = DESeqDataSetFromMatrix(
  countData = countdata,
  colData = sampleInfo,
  design = ~ condition)
dds



countdata <- as.matrix(countdata)
head(countdata)

# Assign condition 
condition <- dds$condition

# Analysis with DESeq2 ----------------------------------------------------

library(DESeq2)

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)

# Plot dispersions
png("qc-dispersions.png", 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Regularized log transformation for clustering/heatmaps, etc
rld <- rlogTransformation(dds)
head(assay(rld))
png("histogram.png", 1000, 1000, pointsize=20)
hist(assay(rld))

dev.off()
# Colors for plots below
## Ugly:
## (mycols <- 1:length(unique(condition)))
## Use RColorBrewer, better
library(RColorBrewer)
(mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))])

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
png("qc-heatmap-samples.png", w=1000, h=1000, pointsize=20)
heatmap.2(as.matrix(sampleDists), key=F, trace="none",
          col=colorpanel(100, "black", "white"),
          ColSideColors=mycols[condition], RowSideColors=mycols[condition],
          margin=c(10, 10), main="Sample Distance Matrix")
dev.off()

# Principal components analysis
## Could do with built-in DESeq2 function:
## DESeq2::plotPCA(rld, intgroup="condition")
## I like mine better:
rld_pca <- function (rld, intgroup = "condition", ntop = 500, colors=NULL, legendpos="bottomleft", main="PCA Biplot", textcx=1, ...) {
  require(genefilter)
  require(calibrate)
  require(RColorBrewer)
  rv = rowVars(assay(rld))
  select = order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(rld)[select, ]))
  fac = factor(apply(as.data.frame(colData(rld)[, intgroup, drop = FALSE]), 1, paste, collapse = " : "))
  if (is.null(colors)) {
    if (nlevels(fac) >= 3) {
      colors = brewer.pal(nlevels(fac), "Paired")
    }   else {
      colors = c("black", "red")
    }
  }
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC1 (",as.character(pc2var),"%)")
  plot(PC2~PC1, data=as.data.frame(pca$x), bg=colors[fac], pch=21, xlab=pc1lab, ylab=pc2lab, main=main, ...)
  with(as.data.frame(pca$x), textxy(PC1, PC2, labs=rownames(as.data.frame(pca$x)), cex=textcx))
  legend(legendpos, legend=levels(fac), col=colors, pch=20)
  #     rldyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$rld),
  #            pch = 16, cerld = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours),
  #                                                                                         terldt = list(levels(fac)), rep = FALSE)))
}
png("qc-pca.png", 1000, 1000, pointsize=20)
rld_pca(rld, colors=mycols, intgroup="condition", xlim=c(-75, 35))
dev.off()


# Get differential expression results
res <- results(dds)
dds
table(res$padj<0.05)
## Order by adjusted p-value
res <- res[order(res$padj), ]
## Merge with normalized count data
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

## Write results
write.csv(resdata, file="diffexpr-results.csv")



## Examine plot of p-values
hist(res$pvalue, breaks=50, col="grey")

## Examine independent filtering
attr(res, "filterThreshold")
plot(attr(res,"filterNumRej"), type="b", xlab="quantiles of baseMean", ylab="number of rejections")

## MA plot
## Could do with built-in DESeq2 function:
## DESeq2::plotMA(dds, ylim=c(-1,1), cex=1)
## I like mine better:
maplot <- function (res, thresh=0.05, labelsig=TRUE, textcx=1, ...) {
  with(res, plot(baseMean, log2FoldChange, pch=20, cex=.5, log="x", ...))
  with(subset(res, padj<thresh), points(baseMean, log2FoldChange, col="red", pch=20, cex=1.5))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<thresh), textxy(baseMean, log2FoldChange, labs=Gene, cex=textcx, col=2))
  }
}
png("diffexpr-maplot.png", 1500, 1000, pointsize=20)
maplot(resdata, main="MA Plot")
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png("diffexpr-volcanoplot.png", 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()

##write results into tables
fdr =0.05
de_total = resdata[which(resdata$padj < fdr),]
length(de_total[,1])
head(de_total)
#Merge the data with biomart file

biomart <-read.csv("Mart_transcripts_genesgr38.p6.94.csv")
de_total <- merge(de_total,biomart,by="Gene")
head(de_total)
write.table((de_total[order(de_total$padj, decreasing=F),]),"DESeq2.de_all.tsv", sep="\t", quote=FALSE, col.names=NA)

de_up = de_total[which(de_total$log2FoldChange>0),]
de_up = de_up[order(de_up$padj, decreasing=F),]   # order by adjusted p-value
length(de_up[,1])
write.table((de_up),"DESeq2.de_up.tsv", sep="\t", quote=FALSE, col.names=NA)

de_down = de_total[which(de_total$log2FoldChange<0),]
de_down = de_down[order(de_down$padj, decreasing=F),]           # order by adjusted p-value
length(de_down[,1])
write.table((de_down),"DESeq2.de_down.tsv", sep="\t", quote=FALSE, col.names=NA)


#new kind of MA plot

jpeg(
  "DESeq2_MAplot.jpeg",
  width=8,
  height=8,
  units="in",
  res=500)
DESeq2::plotMA(dds, alpha=fdr, ylim=c(-2,2),main=sprintf("MA-plot\n(FDR: %.2f, up: %d, down: %d)",fdr,length(de_up[,1]),length(de_down[,1])),ylab="log2 fold change")
dev.off()

# topN genes by pvalue

topN=50
if (length(de_total[,2]) > 0) {
  d = data.frame(id=de_total$GeneName, padj=de_total$padj )
  if ( length(rownames(d)) < topN ) topN = length(rownames(d))
  d
  d_topx_padj = d[order(d$padj, decreasing=F),][1:topN,]
  d_topx_padj
  plotdata = assay(rld)[as.character(d_topx_padj$id),]  # <- error
  plotdata
  
  
  
  if ( exists("gene_names_dic") ) rownames(plotdata) = id_to_gene_name(rownames(plotdata))  # exchange ids by gene names
  plotdata
  
  pdf(sprintf("Fig7.gene_clustering_top%i_DE_genes.pdf",topN), pointsize = 9)
  heatmap.2(plotdata, scale="row", trace="none", dendrogram="column",
            col=colorRampPalette(rev(brewer.pal(9,"RdBu")))(255),
            main=sprintf("Top %d DE genes (by p-value)", topN), keysize=1,
            margins = c(8,10),
            cexRow=0.7, cexCol=0.9)
  dev.off()
}

