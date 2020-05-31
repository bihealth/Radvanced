counts <- as.matrix(read.table("Data/covid_data/counts.tsv"))
annot  <- read.csv("Data/covid_data/annotation.all.csv")
covar  <- read.table("Data/covid_data/covariate_file.txt", sep="\t", header=T)
ds <- DESeqDataSetFromMatrix(countData=counts,
  colData=covar, design=~ 0 + group)
ds <- DESeq(ds)
resultsNames(ds)
## [1] "groupA549.IAV"      "groupA549.IAV.mock" "groupA549.mock"     "groupA549.RSV"      "groupA549.RSV.mock" "groupA549.SC2V"     "groupNHBE.mock"    
## [8] "groupNHBE.SC2V"
contrasts <- list(A549.RSV=c("group", "A549.RSV", "A549.RSV.mock"),
                  A549.IAV=c("group", "A549.IAV", "A549.IAV.mock"),
A549.SC2=c("group", "A549.SC2V", "A549.mock"),
                  NHBE.SC2=c("group", "NHBE.SC2V", "NHBE.mock"))
res <- lapply(contrasts, function(cc) { 
   res <- as.data.frame(results(ds, contrast=cc))
   merge(annot, res, by.x="PrimaryID", by.y=0)
library(DESEq2)
library(DESeq2)
ds <- DESeqDataSetFromMatrix(countData=counts,
  colData=covar, design=~ 0 + group)
ds <- DESeq(ds)
resultsNames(ds)
## [1] "groupA549.IAV"      "groupA549.IAV.mock" "groupA549.mock"     "groupA549.RSV"      "groupA549.RSV.mock" "groupA549.SC2V"     "groupNHBE.mock"    
## [8] "groupNHBE.SC2V"
contrasts <- list(A549.RSV=c("group", "A549.RSV", "A549.RSV.mock"),
                  A549.IAV=c("group", "A549.IAV", "A549.IAV.mock"),
A549.SC2=c("group", "A549.SC2V", "A549.mock"),
                  NHBE.SC2=c("group", "NHBE.SC2V", "NHBE.mock"))
res <- lapply(contrasts, function(cc) { 
   res <- as.data.frame(results(ds, contrast=cc))
   merge(annot, res, by.x="PrimaryID", by.y=0)
 })
head(res[[1]])
length(res)
head(res[[1]])
res[[1]][ order(res[[1]]$pvalue), ]
head(res[[1]][ order(res[[1]]$pvalue), ])
head(res$A549.IAV[ order(res$A549.IAV$pvalue), ])
all(res$A549.SC2$PrimaryID == res$NHBE.SC2$PrimaryID)
sign.a549 <- res$A549.SC2$padj < .01 & abs(res$A549.SC2$log2FoldChange) > 1
sum(sign.a549)
sum(sign.a549, na.rm=T)
head(sign.a549)
sign.a549 <- res$A549.SC2$padj < .05 & abs(res$A549.SC2$log2FoldChange) > 1
sum(sign.a549, na.rm=T)
sign.a549 <- res$A549.SC2$padj < .01 & abs(res$A549.SC2$log2FoldChange) > .5
sum(sign.a549, na.rm=T)
sing.nhbe <- with(res$NHBE.SC2, padj < .01 & abs(log2FoldChange) > .5)
sum(sing.nhbe, na.rm=T)
sum(sing.nhbe & sign.a549)
sum(sing.nhbe & sign.a549, na.rm=T)
plot(res$A549.SC2$log2FoldChange, res$NHBE.SC2$log2FoldChange, pch=19, col="grey", bty="n")
sel <- sing.nhbe & sign.a549
points(res$A549.SC2$log2FoldChange[sel], res$NHBE.SC2$log2FoldChange[sel], pch=19, col="red")
plot(res$A549.SC2$log2FoldChange, res$A549.RSV$log2FoldChange, pch=19, col="grey", bty="n")
sign.a549.rsv <- with(res$A549.RSV, padj < .01 & abs(log2FoldChange) > .5)
sum(sign.a549.rsv, na.rm=TRUE)
sel <- sign.a549.rsv & sign.a549
points(res$A549.SC2$log2FoldChange[sel], res$A549.RSV$log2FoldChange[sel], pch=19, col="red")
library(tmod)
l <- res$A549.SC2$SYMBOL[ order(res$A549.SC2$pvalue) ])
l <- res$A549.SC2$SYMBOL[ order(res$A549.SC2$pvalue) ]
head(l)
l <- res$A549.SC2[ order(res$A549.SC2$pvalue),  ]
head(res$A549.SC2[ order(res$A549.SC2$pvalue),  ])
l <- res$A549.SC2$SYMBOL[ order(res$A549.SC2$pvalue) ]
head(l)
tmodCERNOtest(l)
res_tmod <- tmodCERNOtest(l)
head(res_tmod)
data(tmod)
head(tmod$GENES)
dim(tmod$GENES)
c("1", "23983", "2")
factor(c("1", "23983", "2"))
as.numeric(factor(c("1", "23983", "2")))
getModuleMembers("LI.M75")
getModuleMembers(c("LI.M127", "LI.M75"))
modOverlaps(c("LI.M127", "LI.M75"))
modOverlaps(c("LI.M127", "LI.M75"), stat="n")
intersect(getModuleMembers("LI.M127")[[1]], getModuleMembers("LI.M75")[[1]])
length(intersect(getModuleMembers("LI.M127")[[1]], getModuleMembers("LI.M75")[[1]]))
evidencePlot(l, "LI.M127")
evidencePlot(l, "LI.M75")
evidencePlot(l, "LI.M75", gene.labels=T)
col.v <- rep("grey", length(l))
col.v(sign.a549) <- "red"
col.v[sign.a549] <- "red"
head(col.v)
evidencePlot(l, "LI.M75", gene.labels=T, gene.colors=col.v)
names(col.v) <- l
evidencePlot(l, "LI.M75", gene.labels=T, gene.colors=col.v)
head(col.v)
sum(col.v == "red")
res$A549.SC2[ res$A549.SC2$SYMBOL %in% getModuleMembers("LI.M75")[[1]], ]
col.v[ "OAS1" ]
which(l == "OAS1")
l[8]
col.v[8]
sign.a549[8]
sign.a549 <- res$A549.SC2$padj < .01 & abs(res$A549.SC2$log2FoldChange) > .5
sign.a549[8]
res$A549.SC2[ 8, ]
names(col.v) <- res$A549.SC2$SYMBOL
col.v["OAS1"]
evidencePlot(l, "LI.M75", gene.labels=T, gene.colors=col.v)
res_tmod_full <- lapply(res, function(.r) tmodCERNOtest(.r$SYMBOL[ .r$pvalue ]))
head(res[[1]])
.r <- res[[1]]
head(tmodCERNOtest(.r$SYMBOL[ .r$pvalue ]))
head(.r$SYMBOL[ .r$pvalue ])
class(res)
tmodCERNOtest(.r$SYMBOL[ .r$pvalue ])
l <- .r$SYMBOL[ .r$pvalue ]
tmodCERNOtest(l)
l <- res$A549.SC2$SYMBOL[ order(res$A549.SC2$pvalue) ]
tmodCERNOtest(l)
.r <- res[[1]]
.r$SYMBOL[ order(.r$pvalue) ]
l <- .r$SYMBOL[ order(.r$pvalue) ]
tmodCERNOtest(l)
.r <- res[[1]]
l <- .r$SYMBOL[ order(.r$pvalue) ]
tmodCERNOtest(l)
res_tmod_full <- lapply(res, function(.r) tmodCERNOtest(.r$SYMBOL[ .r$pvalue ]))
res_tmod_full <- lapply(res, function(.r) tmodCERNOtest(.r$SYMBOL[ order(.r$pvalue) ]))
head(res_tmod_full[[1]])
head(res_tmod_full[[2]])
head(res_tmod_full[[3]])
head(res_tmod_full[[3]])
tmodPanelPlot(res)
head(res$A549.RSV)
tmodPanelPlot(res_tmod_full)
tmodPanelPlot(res_tmod_full, filter.rows.pval= 1e-3)
tmodPanelPlot(res_tmod_full, filter.rows.pval= 1e-3, filter.rows.auc=.65)
tmodPanelPlot(res_tmod_full, filter.rows.pval= 1e-3, filter.rows.auc=.8)
lfcs <- sapply(res, function(x) x$log2FoldChange)
head(lfcs)
pvals <- sapply(res, function(x) x$padj)
head(pvals)
lfcs[ is.na(lfcs) ] <- 0
pvals[ is.na(pval) ] <- 1
pvals[ is.na(pvals) ] <- 1
head(pvals)
pie <- tmodDecideTests(res$A549.SC2$SYMBOL, lfc=lfcs, pval=pvals)
head(pie)
head(pie$A549.IAV)
pie <- tmodDecideTests(res$A549.SC2$SYMBOL, lfc=lfcs, pval=pvals, pval.thr=.01, lfc.thr=.5)
tmodPanelPlot(res_tmod_full, filter.rows.pval= 1e-3, filter.rows.auc=.8, pie=pie, grid="b")
