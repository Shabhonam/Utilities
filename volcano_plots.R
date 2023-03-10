library(ggplot2)


volcanoplot <- function(res, lfcthresh=1, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
    with(res, plot(logFC, -log10(PValue), pch=20, main=main, ...))
    with(subset(res, PValue<sigthresh ), points(logFC, -log10(PValue), pch=20, col="red", ...))
    with(subset(res, abs(logFC)>lfcthresh), points(logFC, -log10(PValue), pch=20, col="orange", ...))
    with(subset(res, PValue<sigthresh & abs(logFC)>lfcthresh), points(logFC, -log10(PValue), pch=20, col="green", ...))
    if (labelsig) {
        require(calibrate)
        with(subset(res, PValue<sigthresh & abs(logFC)>lfcthresh), textxy(logFC, -log10(PValue), labs=Gene, cex=textcx, ...))
    }
    legend(legendpos, xjust=1, yjust=1, legend=c(paste("PValue<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}


volcano_plots <- function(CP0, contrast, pval, outdir) {
    common <- intersect(rownames(contrast), CP0$feature_name)
    res_common <- contrast[rownames(contrast) %in% common,]

    png(paste(outdir,"/diffexpr-volcanoplot_",pval,".png",sep=""), 1200, 1000, pointsize=20)
    volcanoplot(contrast, lfcthresh=0.5, sigthresh=pval, textcx=.8, xlim=c(-2.3, 2))
    abline(v=0.5,col="gray",lty=2)
    abline(v=-0.5,col="gray",lty=2)
    abline(v=-log2(1),col="gray",lty=2)
    abline(v=log2(1),col="gray",lty=2)
    abline(h=-log10(pval),col="gray",lty=2)
    dev.off()

    if (length(common)!=0) {
    png(paste(outdir,"/diffexpr-volcanoplot_CP0overlapping_",pval,".png",sep=""), 1200, 1000, pointsize=20)
    volcanoplot(res_common, lfcthresh=0.5, sigthresh=pval, textcx=.8, xlim=c(-2.3, 2))
    abline(v=0.5,col="gray",lty=2)
    abline(v=-0.5,col="gray",lty=2)
    abline(v=-log2(1),col="gray",lty=2)
    abline(v=log2(1),col="gray",lty=2)
    abline(h=-log10(pval),col="gray",lty=2)
    dev.off()
    }
}

directories <- c('GENEHANCERS-Buffy_CoatVSCRC_cfDNA','GENEHANCERS-CRC_tissueVSBuffy_Coat','GENEHANCERS-CRC_tissueVSNAT_tissue','GENES-Buffy_CoatVSCRC_cfDNA','GENES-CRC_tissueVSBuffy_Coat','GENES-CRC_tissueVSNAT_tissue')

#directories <- c("GENEHANCERS-CRC_tissueVSNAT_tissue_2")

for (i in directories) {
    print(i)
    df <- read.csv(paste(i,"results.tsv",sep='/'),sep='\t')
    if (grepl("GENEHANCERS",i,fixed=TRUE)) {
	cp0=read.csv("~/ClassProblem0_top_ml_GENEHANCER_foldC_orig.txt",sep='\t')
	} else {cp0=read.csv("~/ClassProblem0_sig_ml_GENE_foldC.txt",sep='\t')}
    volcano_plots(cp0,df,0.05,i)
}
