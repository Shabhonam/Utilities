library(edgeR)
library(methods)


args <- commandArgs(trailingOnly=TRUE)
case <- args[1]
control <- args[2]
directory <- args[3]
datasource <- args[4]
features <- args[5]

df <- read.csv(datasource,sep='\t')
feat <- read.delim(features,row.names=1,sep='\t')
PooledCounts <- sumTechReps(feat, ID=df$to_collapse)
df <- df[!duplicated(df$to_collapse), ]

if (all(order(c(control,case))==c(1,2))) {
    df <- df[order(df$patient_id,df$comparisons),]
} else {
    df <- df[order(df$patient_id,rev(df$comparisons)),]
}

dge <- DGEList(counts=PooledCounts,samples=df)
design <- model.matrix(~dge$samples$patient_id+dge$samples$comparisons)
disp <- estimateDisp(dge,design)

fit <- glmQLFit(disp, design)
qlf <- glmQLFTest(fit)
write.table(topTags(qlf,n="Inf"),file=paste(directory,"/results.tsv",sep=''),sep='\t')
