#make sure to load R/3.6.2

library( "DESeq2" )
library( "gplots" )
library( "RColorBrewer" )
library( "genefilter" )
library( "EnhancedVolcano" )

sampleInfo <- read.table("../shortRNAseq.txt")
sampleInfo$FullSampleName <- as.character(sampleInfo$FullSampleName)

countdata <- read.table("fly_counts_2.txt", header=TRUE, row.names=1)
countdata <- countdata[ ,6:ncol(countdata)]

temp <- colnames(countdata)
temp <- gsub("RNAseq.bam.","",temp)
temp <- gsub(".sorted.bam","",temp)
temp <- gsub("...alignments.","",temp)
colnames(countdata) <- temp

sampleInfo <- sampleInfo[order(sampleInfo$FullSampleName),] 
cbind(temp,sampleInfo$FullSampleName,temp == sampleInfo$FullSampleName)

# create DEseq2 object & run DEseq
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=sampleInfo, design=~TissueCode)
dds <- DESeq(dds)
res <- results( dds )

plotMA( res, ylim = c(-1, 1) )
plotDispEsts( dds )
hist( res$pvalue, breaks=20, col="grey" )
res <- res[res$baseMean>0,] 
rld <- rlog( dds )

sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- rld$TissueCode
colnames(sampleDistMatrix) <- NULL
colours <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
print( plotPCA( rld, intgroup = "TissueCode") )
topVarGenes <- head( order( rowVars( assay(rld) ), decreasing=TRUE ), 35 )
heatmap.2( assay(rld)[ topVarGenes, ], scale="row", trace="none", dendrogram="column", col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

# volcano plot
#from https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html

res <- results(dds,
    contrast = c('TissueCode','E','B'))
res <- lfcShrink(dds,
    contrast = c('TissueCode','E','B'), res=res, type = 'normal')
res <- res[res$baseMean>0,]

EnhancedVolcano(res,
    lab = rownames(res),
    x = 'log2FoldChange',
    y = 'pvalue')
