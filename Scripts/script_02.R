##------ Wed May  6 10:18:58 2020 ------##
foo <- read.table("../Radvanced/Data/covid_data/counts.tsv")
head(foo)
plot(density(foo[,1]))
fool <- load("Data/covid_data/DESeq2.all.rld.model.csv")
fool <- read.table("Data/covid_data/DESeq2.all.rld.model.csv")
head(fool)
fool <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE)
head(fool)
plot(foo[,1], foo[,2], log="x")
plot(foo[,1], foo[,2])
library(GEOquery)
?getGEO
?getGEO
foo <- getGEO("GSE147507")
foo
class(foo[[1]])
foo[[1]]
exprs(foo[[1]])
exprs(foo[[2]])
ls()
os()
foo[[2]]
counts <- read.table("Data/covid_data/counts.tsv")
head(counts)
annotation <- read.csv("Data/covid_data/annotation.all.csv")
head(annotation)
all(rownames(counts) == annotation$PrimaryID)
head(annotation$SYMBOL[ match(rownames(counts), annotation$PrimaryID)])
foo <- read.table("Data/covid_data/contrast.A549_IAV_ID3.csv")
head(foo)
foo <- read.csv("Data/covid_data/contrast.A549_IAV_ID3.csv")
head(foo)
foo <- read.table("Data/covid_data/contrast.A549_IAV_ID3.csv", header=TRUE)
head(foo)
head(foo[ order(foo$pvalue), ])
colnames(counts)
cell.types <- c(rep("NHBE", 6), rep("A549", 14))
cell.types
treatment <- c(rep(c("M", "S2"), each=3), rep(c("M", "S2"), each=3), rep(c("M", "RSV"), each=2), rep(c("M", "IAV"), each=2))
length(treatment)
length(cell.types)
covariate_table <- data.frame(ID=colnames(counts), cell.type=cell.types, treatment=treatment)
covariate_table
head(foo)
head(foo[ order(foo$pvalue), ])
foo$Description <- annotation$GENENAME[ match(foo$gene_id, annotation$PrimaryID)]
head(foo)
head(foo[ order(foo$pvalue), ])
head(foo[ order(foo$pvalue), ], 20)
foo <- read.table("Data/covid_data/contrast.A549_RSV_ID2.csv", header=TRUE)
foo$Description <- annotation$GENENAME[ match(foo$gene_id, annotation$PrimaryID)]
head(foo[ order(foo$pvalue), ], 20)
foo.m <- merge(foo, annotation, by.x="gene_id", by.y="PrimaryID")
head(foo.m)
head(counts)
mtx <- data.matrix(counts)
head(mtx)
class(mtx)
head(foo)
head(foo[ order(foo$pvalue), ], 20)
id <- "ENSG00000115415"
boxplot(counts[id, ] ~ covariate_table$cell.type)
boxplot(mtx[id, ] ~ covariate_table$cell.type)
counts[id,]
class(counts[id,])
mtx <- data.matrix(counts)
class(mtx[id,])
mtx[id,]
boxplot(as.vector(counts[id, ]) ~ covariate_table$cell.type)
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type)
plot(density(counts[,1]))
plot(density(log2(counts[,1])))
lines(density(log2(counts[,2])), col="red")
lines(density(log2(counts[,3])), col="blue")
lines(density(log2(counts[,4])), col="green")
d <- apply(counts, 2, density)
length(d)
class(d)
class(d[[1]])
plot(d[[1]])
d <- apply(log2(counts), 2, density)
plot(NULL, xlim=c(0, max(log2(counts))))
plot(NULL, xlim=c(0, max(log2(counts))), ylim=c(0, 0.1))
lines(d[[1]])
plot(NULL, xlim=c(0, max(log2(counts))), ylim=c(0, 0.1))
lapply(1:20, function(i) lines(d[[i]], col=i))
plot(NULL, xlim=c(0, max(log2(counts))), ylim=c(0, 0.2))
lapply(1:20, function(i) lines(d[[i]], col=i))
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type)
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type, log="y")
library(beeswarm)
beeswarm(unlist(counts[id, ]) ~ covariate_table$cell.type, log="y")
beeswarm(unlist(counts[id, ]) ~ covariate_table$cell.type)
beeswarm(unlist(counts[id, ]) ~ covariate_table$cell.type, log=T)
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type, log="y", add=T)
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type, log="y")
beeswarm(unlist(counts[id, ]) ~ covariate_table$cell.type, log=T, add=T)
head(foo[ order(foo$pvalue), ], 20)
boxplot(unlist(counts[id, ]) ~ covariate_table$treatment, log="y")
covariate_table$c.t <- paste0(covariate_table$cell.type, ".", covariate_table$treatment)
covariate_table
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y")
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y", text.cex=.7)
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y", cex=.7)
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y", las=90)
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y", las=1)
boxplot(unlist(counts[id, ]) ~ covariate_table$c.t, log="y", las=2)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T)
counts2 <- read.table("DESeq2.all.rld.model.csv")
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv")
head(counts2)
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE)
head(counts2)
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1)
head(counts2)
boxplot(unlist(counts2[id, ]) ~ covariate_table$c.t, log="y", las=2)
colnames(counts2)
counts2 <- counts2[ , 1:20]
boxplot(unlist(counts2[id, ]) ~ covariate_table$c.t, log="y", las=2)
colSums(ccounts)
colSums(counts2)
colSums(counts)
colSums(counts) / colSums(counts2)
colSums(counts2) / colSums(counts)
beeswarm(unlist(counts2[id, ]) ~ covariate_table$c.t, add=TRUE)
boxplot(unlist(counts[id, ]) ~ covariate_table$treatment, LOG
asdjfad
´
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1)
counts2 <- read.csv("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1, sep="\t")
head(counts2)
counts2 <- read.csv("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, sep="\t")
head(counts2)
counts2 <- read.csv("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, sep="\t", row.names=1)
counts2 <- read.csv("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, sep="\t", row.names=1)[,1:20]
annotation[ grep("IFNG", annotation$GENENAME), ]
annotation[ grep("IFN", annotation$GENENAME), ]
annotation[ grep("IFNG", annotation$SYMBOL), ]
id <- "ENSG00000111537"
boxplot(unlist(counts2[id, ]) ~ covariate_table$treatment, las=2)
annotation[ grep("CXCL5", annotation$SYMBOL), ]
id <- "ENSG00000163735"
boxplot(unlist(counts2[id, ]) ~ covariate_table$treatment, las=2)
beeswarm(unlist(counts2[id, ]) ~ covariate_table$c.t, add=TRUE)
boxplot(unlist(counts2[id, ]) ~ covariate_table$c.t, las=2)
beeswarm(unlist(counts2[id, ]) ~ covariate_table$c.t, add=TRUE)
##------ Fri May  8 13:11:15 2020 ------##
counts
d <- apply(log2(counts), 2, density)
plot(NULL, xlim=c(0, max(log2(counts))), ylim=c(0, 0.2))
lapply(1:20, function(i) lines(d[[i]], col=i))
d[[1]]
d[[1]]$x
d[[2]]$x
d.gg <- lapply(d, function(x) cbind(x$x, x$y))
d.gg[[1]]
d.gg <- lapply(1:ncol(counts), function(i) cbind(d[[i]]$x, d[[i]]x$y, colnames(counts)[i]))
d.gg <- lapply(1:ncol(counts), function(i) cbind(d[[i]]$x, d[[i]]$y, colnames(counts)[i]))
head(d.gg[[1]])
d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=d[[i]]$x, y=d[[i]]x$y, ID=colnames(counts)[i]))
d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=d[[i]]$x, y=d[[i]]$y, ID=colnames(counts)[i]))
head(d.gg[[1]])
d.gg <- Reduce(rbind, d.gg)
dim(d.gg)
head(d.gg)
ggplot(d.gg, aes(x=x, y=y)) + geom_l
library(ggplot2)
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line()
remove.packages("farver")
install.packages("farver")
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line()
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line(col=group)
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line(aes(col=group))
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line(aes(col=ID))
d.gg <- lapply(1:ncol(counts), function(i) data.frame(y=counts[,i], ID=colnames(counts)[i]))
d.gg <- Reduce(rbind, d.gg)
g
dim8d.gg)
dim(d.gg)
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_density()
ggplot(d.gg, aes(y=y, group=ID)) + geom_density()
d.gg$y <- log2(d.gg$y)
ggplot(d.gg, aes(y=y, group=ID)) + geom_density()
head(d.gg)
ggplot(d.gg, aes(x=y, group=ID)) + geom_density()
ggplot(d.gg, aes(x=y, group=ID)) + geom_density(aes(col=ID))
d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=counts[,i], ID=colnames(counts)[i]))
d.gg <- Reduce(rbind, d.gg)
d.gg$y <- log2(d.gg$y)
## check how large this data frame is!
dim(d.gg)
ggplot(d.gg, aes(x=x, group=ID)) + geom_density(aes(col=ID))
d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=counts[,i], ID=colnames(counts)[i]))
d.gg <- Reduce(rbind, d.gg)
d.gg$x <- log2(d.gg$x)
## check how large this data frame is!
dim(d.gg)
ggplot(d.gg, aes(x=x, group=ID)) + geom_density(aes(col=ID))
gg <- data.frame(covariate_table, x=mtx[id,])
gg
mtx <- data.matrix(counts)
## let us take a look at one of the top genes in the contrast
id <- "ENSG00000115415"
## this will not work: counts is a data frame and so counts[id, ] is also a
## data frame.
boxplot(counts[id, ] ~ covariate_table$cell.type)
## we can either do this
boxplot(unlist(counts[id, ]) ~ covariate_table$cell.type)
## or this
boxplot(mtx[id, ] ~ covariate_table$cell.type)
## but actually, we should log the data!
boxplot(mtx[id, ] ~ covariate_table$cell.type, log="y")
## we can add a beeswarm
library(beeswarm)
beeswarm(unlist(counts[id, ]) ~ covariate_table$cell.type, log=T, add=T)
## but actually, we should log the data!
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y")
## we can add a beeswarm
library(beeswarm)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", bty="n")
?boxplot
?boxplot
?bxp
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", bty="n", frame.plot="n")
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", bty="n", frame.plot=F)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, main=id)
## make a nicer version
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, ylab="Raw counts", lab.cex=.7)
  main=id, ylab="Raw counts", lab.cex=.2)
## make a nicer version
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, ylab="Raw counts", lab.cex=.7)
## make a nicer version
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, ylab="Raw counts", lab.cex=.5)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, lab.cex=2)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, cex.lab=2)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, cex.axis=2)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, cex.axis=.5)
library(RColorBrewer)
display.brewer.all()
brewer.pal("dark2")
brewer.pal(7, "dark2")
brewer.pal(7, "Dark2")
length(unique(covariate_table$c.t))
xclip(15)
covariate_table$color <- brewer.pal(6, "Dark2")[ factor(covariate_table$c.t) ]
covariate_table
?beeswarm
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE,
  main=id, ylab="Raw counts", cex.axis=.7)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T,
  col=brewer.pal(6, "Dark2"))
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE,
  main=id, ylab="Raw counts", cex.axis=.7)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T,
  col=brewer.pal(6, "Dark2"), pch=19)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T,
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)
ggplot(gg, aes(x=x, group=c.t)) + geom_boxplot()
ggplot(gg, aes(x=x, group=c.t)) + geom_boxplot(aes(col=c.t))
ggplot(gg, aes(x=x, group=c.t)) + geom_boxplot(aes(col=c.t)) + scale_x_continuous(trans=log2_trans())
library(scales)
ggplot(gg, aes(x=x, group=c.t)) + geom_boxplot(aes(col=c.t)) + scale_x_continuous(trans=log2_trans())
ggplot(gg, aes(x=c.t, y=x)) + geom_boxplot(aes(col=c.t)) + scale_x_continuous(trans=log2_trans())
ggplot(gg, aes(x=c.t, y=x)) + geom_boxplot(aes(col=c.t)) + scale_y_continuous(trans=log2_trans())
xclip(1)
gg <- data.frame(covariate_table, y=mtx[id,])
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + scale_y_continuous(trans=log2_trans())
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm() + scale_y_continuous(trans=log2_trans()) 
??geom_beeswwarm
??geom_beeswarm
install.packages("ggbeeswarm")
library("ggbeeswarm")
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm() + scale_y_continuous(trans=log2_trans()) 
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t)) + scale_y_continuous(trans=log2_trans()) 
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t), cex=1.5) + scale_y_continuous(trans=log2_trans()) 
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t), size=1.5) + scale_y_continuous(trans=log2_trans()) 
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t), size=3.5) + scale_y_continuous(trans=log2_trans()) 
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t), size=2.5) + scale_y_continuous(trans=log2_trans()) 
library(cowplot)
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) + geom_beeswarm(aes(col=c.t), size=2.5) + scale_y_continuous(trans=log2_trans()) 
?pairs
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1)
mtx2 <- data.matrix(counts2)
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1)
mtx2 <- data.matrix(counts2[, 1:20])
head(mtx2)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, ylab="Raw counts", cex.axis=.7)
beeswarm(mtx[id, ] ~ covariate_table$c.t, log=T, add=T, 
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)
boxplot(mtx2[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE,
  main=id, ylab="Raw counts", cex.axis=.7)
beeswarm(mtx2[id, ] ~ covariate_table$c.t, log=T, add=T,     
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)
boxplot(mtx2[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, xlab="Experimental group", ylab="Raw counts", cex.axis=.7)
beeswarm(mtx2[id, ] ~ covariate_table$c.t, log=T, add=T,     
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)
##------ Tue May 12 10:51:38 2020 ------##
list.files("Data/covid_data/", pattern="contrast.*")
list.files("Data/covid_data/", pattern="contrast*")
read.csv("Data/covid_data/contrast.A549_IAV_ID3.csv")
head(read.table("Data/covid_data/contrast.A549_IAV_ID3.csv"))
head(read.table("Data/covid_data/contrast.A549_IAV_ID3.csv", header=TRUE))
    files.list <- list.files("Data/covid_data/", pattern="contrast*")
    files <- lapply(files.list, read.table, header=TRUE, )
files.list
?list.files
    files.list <- list.files("Data/covid_data/", pattern="contrast*", full.names=TRUE)
files.list
    files.list <- list.files("Data/covid_data/", pattern="contrast*", full.names=TRUE)
    files <- lapply(files.list, read.table, header=TRUE)
files
    data <- lapply(files, function(x) x$log2FoldChange)
    # or
    data <- lapply(files, function(x) x[["log2FoldChange"]])
head(data[[1]])
    data <- lapply(files, `[[`, "log2FoldChange")
head(data[[1]])
    mtx <- Reduce(cbind, data)
head(mtx)
    mtx <- sapply(data, `[[`, "log2FoldChange")
    mtx <- sapply(files, `[[`, "log2FoldChange")
head(mtx)
    `[[`(files, "log2FoldChange")
files
    `[[`(files, 1)
    plot(mtx[,1], mtx[,2], pch=19)
files
lapply(files, nrow)
head(files)
head(files[[1]])
head(files[[2]])
    genes <- files[[1]]$gene_id
    files <- lapply(files, function(x) x[ match(genes, x$gene_id), ])
head(files[[2]])
head(files[[3]])
    mtx <- sapply(files, `[[`, "log2FoldChange")
class(mtx)
dim(mtx)
    plot(mtx[,1], mtx[,2], pch=19)
    plot(mtx[,1], mtx[,3], pch=19)
    plot(mtx[,1], mtx[,4], pch=19)
install.packages("corrgram")
remove.packages("gtools")
install.packages("gtools")
install.packages("corrgram")
plot(mtx)
plot(as.data.frame(mtx))
plot(as.data.frame(mtx), pch=19, col="#33333333")
    cor(mtx)
    cor(mtx, na.rm=TRUE)
    cor(mtx, use="p")
?cor
dim(mtx)
dim(files[[1]])
library(corrgram)
?corrgram
corrgram(mtx)
corrgram(mtx, lower.panel=panel.conf)
corrgram(mtx, lower.panel=panel.pie)
corrgram(mtx, lower.panel=panel.vote)
corrgram(mtx, lower.panel=panel.minmax)
corrgram(mtx, lower.panel=panel.pts)
corrgram(mtx, lower.panel=panel.pts, upper.panel=panel.density)
corrgram(mtx, lower.panel=panel.ellipse, upper.panel=panel.cor)
corrgram(mtx, lower.panel=panel.ellipse, upper.panel=panel.shade)
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse)
?corrgram
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.fill)
?corrgram
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.bar)
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse)
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse, diag.panel=panel.txt)
?corrgram
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse, diag.panel=panel.minmax)
corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse, diag.panel=panel.density)
library(ellipse)
?pairs
plotcorr(mtx)
plotcorr(cor(mtx))
?plotcorr
plotcorr(cor(mtx), col="grey")
install.packages("GGAlly")
install.packages("GGally")
library(GGally)
ggpairs(mtx)
ggpairs(data.frame(mtx))
ggpairs(data.frame(mtx), method=c("everything", "pearson"))
?pairs
?pairs
pairs(mtx)
?smoothScatter
    u.panel.f <- function(x, y) {
      smoothScatter(x, y)
    }
    pairs(mtx, upper.panel=u.panel.f)
?smoothScatter
    u.panel.f <- function(x, y) {
      smoothScatter(x, y, add=TRUE)
    }
    pairs(mtx, upper.panel=u.panel.f)
?pairs
foo <- cor(files[,1], files[,2]))
foo <- cor(mtx[,1], mtx[,2])
foo
foo <- cor.test(mtx[,1], mtx[,2], use="p")
foo
foo$statistic
foo$estimate
format.pval(foo$p.value)
format.pval(foo$p.value, 2)
?text
    u.panel.f <- function(x, y) {
      smoothScatter(x, y, add=TRUE)
    }
    l.panel.f <- function(x, y) {
      ctst <- cor.test(x, y, use="p")
      par(usr=c(0, 1, 0, 1))
      txt <- sprintf("r²=%.2f\np=%s",
        ctst$estimate,
        format.pval(ctst$p.value, 2))
      text(.5, .5, txt)
    }
    pairs(mtx, upper.panel=u.panel.f, lower.panel=l.panel.f)
    smoothScatter(mtx[,1], mtx[,4])
    cid <- 1
    sel.contrast <- abs(files[[cid]]$log2FoldChange) > 1
    sel.pval     <- files[[cid]]$padj < .01
    sum(sel.contrast & sel.pval, na.rm=TRUE)
cid <- 2
    sel.contrast <- abs(files[[cid]]$log2FoldChange) > 1
    sel.pval     <- files[[cid]]$padj < .01
    sum(sel.contrast & sel.pval, na.rm=TRUE)
    ## this function takes a data frame
    get.stats <- function(x,  pval.thr=.01, lfc.thr=1) {
      sel.contrast <- abs(x$log2FoldChange) > 1
      sel.pval     <- x$padj < .01
      return(sum(sel.contrast & sel.pval, na.rm=TRUE))
    } 
    ctr.stats <- function(x, pval.thr=.01, lfc.thr=1) {
      sapply(x, get.stats)
    }
ctr.stats(files)
    ctr.stats(files, pval.thr=.001)
    ## this function takes a data frame
    get.stats <- function(x,  pval.thr=.01, lfc.thr=1) {
      sel.contrast <- abs(x$log2FoldChange) > 1
      sel.pval     <- x$padj < .01
      return(sum(sel.contrast & sel.pval, na.rm=TRUE))
    } 
    ctr.stats <- function(x, pval.thr=.01, lfc.thr=1) {
      sapply(x, get.stats, pval.thr=pval.thr, lfc.thr=lfc.thr)
    }
    
    ctr.stats(files)
    ctr.stats(files, pval.thr=.001)
    get.stats <- function(x,  pval.thr=.01, lfc.thr=1) {
      sel.contrast <- abs(x$log2FoldChange) > lfc.thr
      sel.pval     <- x$padj < pval.thr
      return(sum(sel.contrast & sel.pval, na.rm=TRUE))
    } 
    ctr.stats <- function(x, pval.thr=.01, lfc.thr=1) {
      sapply(x, get.stats, pval.thr=pval.thr, lfc.thr=lfc.thr)
    }
    
    ctr.stats(files)
    ctr.stats(files, pval.thr=.001)
    ctr.stats(files, pval.thr=.0001)
    ctr.stats(files, pval.thr=.0001, lfc.thr=2)
plot(density(rnorm(1000))
)
?rbin
?rnbinom
?rbinom
?rbinom
plot(density(rbinom(1000, .5)))
?rbinom
?rbinom
plot(density(rbinom(1000, 500, .5)))
plot(density(rbinom(1000, 50, .5)))
plot(density(rbinom(1000, 10, .5)))
plot(density(rbinom(1000, 10, .1)))
library(RColorBrewer)
display.brewer.all
display.brewer.all()
brewer.pal("Set2", n=10)
brewer.pal("Set2", n=8)
library(RColorBrewer)
pal <- brewer.pal("Set2", n=8)
norm <- rnorm(1000)
binom1 <- rbinom(1000, 10, .5)
binom2 <- rbinom(1000, 10, .3)
binom3 <- rbinom(1000, 10, .1)
plot(density(binom1), col=pal[1])
lines(density(binom2), col=pal[2])
plot(density(binom1), col=pal[1], ylim=c(0, .4))
lines(density(binom2), col=pal[2])
plot(density(binom1), col=pal[1], ylim=c(0, .4), lwd=3)
lines(density(binom2), col=pal[2], lwd=3)
lines(density(binom3), col=pal[3], lwd=3)
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3)
lines(density(binom2), col=pal[2], lwd=3)
lines(density(binom3), col=pal[3], lwd=3)
polygon(density(binom3))
polygon(density(binom3), col=pal[3])
?polygon
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3,
  main="Binomial distribution")
polygon(density(binom1), col=pal[1], border=NA)
pal <- brewer.pal("Set2", n=8)
pal.t <- paste0(pal, '33')
norm <- rnorm(1000)
binom1 <- rbinom(1000, 10, .5)
binom2 <- rbinom(1000, 10, .3)
binom3 <- rbinom(1000, 10, .1)
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3,
  main="Binomial distribution")
polygon(density(binom1), col=pal.t[1], border=NA)
polygon(density(binom2), col=pal.t[2], border=NA)
lines(density(binom2), col=pal[2], lwd=3)
polygon(density(binom3), col=pal.t[3], border=NA)
lines(density(binom3), col=pal[3], lwd=3)
pal <- brewer.pal("Set2", n=8)
pal.t <- paste0(pal, '33')
norm <- rnorm(1000)
binom1 <- rbinom(1000, 10, .5)
binom2 <- rbinom(1000, 10, .3)
binom3 <- rbinom(1000, 10, .1)
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3,
  xlab="Number of heads in 10 throws",
  main="Binomial distribution", bty="n")
polygon(density(binom1), col=pal.t[1], border=NA)
polygon(density(binom2), col=pal.t[2], border=NA)
lines(density(binom2), col=pal[2], lwd=3)
polygon(density(binom3), col=pal.t[3], border=NA)
lines(density(binom3), col=pal[3], lwd=3)
legend("topright", c("p = 0.05", "p = 0.03", "p = 0.01"),
  col=pal.t, lwd=3)
pal <- brewer.pal("Set2", n=8)
pal.t <- paste0(pal, '33')
norm <- rnorm(1000)
binom1 <- rbinom(1000, 10, .5)
binom2 <- rbinom(1000, 10, .3)
binom3 <- rbinom(1000, 10, .1)
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3,
  xlab="Number of heads in 10 throws",
  main="Binomial distribution", bty="n")
polygon(density(binom1), col=pal.t[1], border=NA)
polygon(density(binom2), col=pal.t[2], border=NA)
lines(density(binom2), col=pal[2], lwd=3)
polygon(density(binom3), col=pal.t[3], border=NA)
lines(density(binom3), col=pal[3], lwd=3)
legend("topright", c("p = 0.05", "p = 0.03", "p = 0.01"),
  col=pal, lwd=3)
plot(density(binom1), col=pal[1], ylim=c(0, .7), lwd=3,
  xlab="Number of heads in 10 throws",
  main="Binomial distribution", bty="n")
polygon(density(binom1), col=pal.t[1], border=NA)
polygon(density(binom2), col=pal.t[2], border=NA)
lines(density(binom2), col=pal[2], lwd=3)
polygon(density(binom3), col=pal.t[3], border=NA)
lines(density(binom3), col=pal[3], lwd=3)
legend("topright", c("p = 0.05", "p = 0.03", "p = 0.01"),
  col=pal, lwd=3, bty="n")
hist(binom1)
lines(density(binom1))
?hist
lines(density(binom1), freq=TRUE)
hist(binom1, freq=TRUE)
hist(binom1, freq=TRUE)
hist(binom1, freq=F)
lines(density(binom1))
hist(binom1, freq=F, ylim=c(0, .7))
lines(density(binom1))
hist(binom2, freq=F, ylim=c(0, .7), add=TRUE)
?hist
?hist
hist(binom3, col=pal.t[3], border=pal[3], lwd=3)
hist(binom3, col=pal.t[3], border=pal[3], lwd=3)
hist(binom3, col=pal.t[3], border=pal[3])
hist(binom3, col=pal.t[3], border=pal[3])
hist(binom3, col=pal.t[2], border=pal[2], add=TRUE)
hist(binom3, col=pal.t[1], border=pal[1], add=TRUE)
hist(binom3, col=pal.t[3], border=pal[3])
hist(binom2, col=pal.t[2], border=pal[2], add=TRUE)
hist(binom1, col=pal.t[1], border=pal[1], add=TRUE)
hist(binom3, col=pal.t[3], border=pal[3], ylim=c(0, .7), xlim=c(0, 10))
hist(binom2, col=pal.t[2], border=pal[2], add=TRUE)
hist(binom1, col=pal.t[1], border=pal[1], add=TRUE)
hist(binom3, col=pal.t[3], freq=FALSE,  border=pal[3], ylim=c(0, .7), xlim=c(0, 10))
hist(binom2, col=pal.t[2], freq=FALSE,  border=pal[2], add=TRUE)
hist(binom1, col=pal.t[1], freq=FALSE,  border=pal[1], add=TRUE)
?hist
hist(binom3, col=pal.t[3], freq=FALSE,  border=pal[3], ylim=c(0, .7), xlim=c(0, 10))
hist(binom3, col=pal.t[3], freq=FALSE,  border=pal[3], ylim=c(0, .7), xlim=c(0, 10), breaks=1:10)
hist(binom3, col=pal.t[3], freq=FALSE,  border=pal[3], ylim=c(0, .7), xlim=c(0, 10), breaks=0:10)
hist(binom3, col=pal.t[3], freq=FALSE,  border=pal[3], ylim=c(0, .8), xlim=c(0, 10), breaks=0:10)
hist(binom2, col=pal.t[2], freq=FALSE,  border=pal[2], add=TRUE,
breaks=0:10) 
hist(binom1, col=pal.t[1], freq=FALSE,  border=pal[1], add=TRUE,
breaks=0:10)
norms <- lapply(1:3, function(x) rnorm(1000, sd=x))
norms[[1]]
norms.d <- lapply(norms, d)
norms.d <- lapply(norms, density)
plot(norms.d[[1]])
plot(norms.d[[3]])
plot(norms.d[[1]])
plot(norms.d[[3]])
norms <- lapply(1:3, function(x) rnorm(1000, sd=x))
norms.d <- lapply(norms, density)
plot(NULL, xlim=c(-10, 10), ylim=c(0, .5), bty="n")
lapply(1:3, function(x) lines(norms.d[[x]], col=pal[x], lwd=3))
library(RColorBrewer)
pal <- brewer.pal("Set2", n=8)
pal.t <- paste0(pal, '33')
norms <- lapply(1:3, function(x) rnorm(1000, sd=x))
norms.d <- lapply(norms, density)
plot(NULL, xlim=c(-10, 10), ylim=c(0, .5), bty="n")
lapply(1:3, function(x) {
  polygon(norms.d[[x]], col=pal.t[x], border=NA)
  lines(norms.d[[x]], col=pal[x], lwd=3)
})
x <- norms[[1]]
sqrt(sum(x - mean(x))^2/n)
n <- length(x)
sqrt(sum(x - mean(x))^2/n)
sqrt(sum(x - mean(x))^2/n^2)
sqrt(sum((x - mean(x))^2)/n)
counts <- read.table("Data/covid_data/counts.tsv")
head(counts)
annot <- read.csv("Data/covid_data/annotation.all.csv")
head(annot)
covar <- read.table("Data/covid_data/covariate_file.txt", sep="\t")
head(covar)
covar <- read.table("Data/covid_data/covariate_file.txt", sep="\t", header=T)
xclip(5)
counts <- read.table("Data/covid_data/counts.tsv")
xclip(1)
head(covar)
head(counts)
all(colnames(counts) %in% covar$label)
d <- density(counts[,1])
plot(d)
plot(d, ylim=c(0, .001))
plot(d, ylim=c(0, .001), xlim=c(0, 10000))
plot(d, ylim=c(0, .001), xlim=c(0, 1000))
plot(d, ylim=c(0, .001), xlim=c(0, 10000))
plot(d, ylim=c(0, .0001), xlim=c(0, 10000))
xclip(1)
?rnbinom
rnbinom(10, 100, .5)
?rnbinom
rnb <- lapply(c(.5, .3, .1), function(x) rnbinom(1000, 100, x))
rnb <- lapply(c(.5, .3, .1), function(x) rnbinom(1000, 100, x))
rnb.d <- lapply(rnb, density)
plot(NULL, xlim=c(-10, 10), ylim=c(0, .5), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d[[x]], col=pal.t[x], border=NA)
  lines(rnb.d[[x]], col=pal[x], lwd=3)
})
plot(density(rnb.d[[1]]))
rnb.d[[1]]
plot(rnb.d[[1]])
plot(rnb.d[[2]])
plot(rnb.d[[3]])
rnb <- lapply(c(.9, .5, .3, .1), function(x) rnbinom(1000, 100, x))
rnb.d <- lapply(rnb, density)
plot(NULL, xlim=c(-10, 10), ylim=c(0, .5), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d[[x]], col=pal.t[x], border=NA)
  lines(rnb.d[[x]], col=pal[x], lwd=3)
})
plot(NULL, xlim=c(0, 300), ylim=c(0, .5), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d[[x]], col=pal.t[x], border=NA)
  lines(rnb.d[[x]], col=pal[x], lwd=3)
})
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d[[x]], col=pal.t[x], border=NA)
  lines(rnb.d[[x]], col=pal[x], lwd=3)
})
rnb <- lapply(c(.9, .5, .3), function(x) rnbinom(1000, 100, x)) 
rnb.d <- lapply(rnb, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d[[x]], col=pal.t[x], border=NA)
  lines(rnb.d[[x]], col=pal[x], lwd=3)
})
paste0("p = ", c(.9, .5, .3))
length(paste0("p = ", c(.9, .5, .3)))
legend("topright", paste0("p = ", c(.9, .5, .3)), col=pal, lwd=3, bty="n")
BiocManager::installlibrary(DESeq2)
library(DESeq2)
x <- 1:10
y <- x + rnorm(10) # add noise
lm(y ~ x)
x <- 1:10
y <- x + rnorm(10) # add noise
summary(lm(y ~ x))
x <- factor(rep(c("A", "B"), each=5))
y <- rep(1:2, each=5) + rnorm(10)
x
plot(y ~ x)
anova(lm(y ~ x))
aov(lm(y ~ x))
summary(lm(y ~ x))
plot(y ~ x)y <- rep(1:3, each=5) + rnorm(10)
y <- rep(1:3, each=5) + rnorm(10)
y <- rep(c(1,3), each=5) + rnorm(10)
summary(lm(y ~ x))
aov(lm(y ~ x))
anova(lm(y ~ x))
anova(lm(y ~ 0 + x))
summary(lm(y ~ 0 + x))
DESeqDataSetFromMatrix(countData=counts, colData=covar, design="~ group")
DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ group)
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ group)
ds
DESeq(ds)
class(ds)
ds <- DESeq(ds)
class(ds)
resultsNames(ds)
unique(covar$group)
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ 0 + group)
class(ds)
ds <- DESeq(ds)
resultsNames(ds)
?results
covar
covar$cell <- gsub("\\..*", "", covar$group)
covar
covar$treatment <- gsub(".*\\.", "", covar$group)
covar
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ 0 + cell + treatment)
resultsNames(ds)
ds <- DESeq(ds)
resultsNames(ds)
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ 0 + cell * treatment)
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar, design=~ cell * treatment)
model.matrix(~ 0 + cell * treatment, data=covar)
?DESeqDataSetFromMatrix
ds <- DESeqDataSetFromMatrix(countData=counts, colData=covar)
?DESeq
?rnbinom
lapply(rnb, mean)
mus <- sapply(rnb, mean)
mus
?rnbinom
p <- c(.9, .5, .3)
p * mus
p * mus / (1 - p)
mus <- c(10, 100, 233)
alpha <- c(98.94600, 100.07600,  99.62829)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
plot(density(rnb2[[1]]))
plot(density(rnb2[[2]]))
mus <- c(10, 100, 233)
alpha <- c(98.94600, 100.07600,  99.62829)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=100", mus), col=pal, lwd=3, bty="n")
mus <- c(10, 100, 233)
alpha <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=100", mus), col=pal, lwd=3, bty="n")
counts(ds)
ds <- ds[ rs > 10, ]
rs <- rowSums(counts(ds))
ds <- ds[ rs > 10, ]
sum(rs > 10)
length(rs)
cpm <- counts(ds) / colSums(counts(ds))
colSums(cpm)
cpm <- counts(ds) / rep(colSums(counts(ds)), nrow(counts(ds)))
colSums(cpm)
cpm <- counts(ds) / rep(colSums(counts(ds)), each=nrow(counts(ds)))
colSums(cpm)
cpm <- counts(ds) / rep(colSums(counts(ds)), each=nrow(counts(ds))) * 1e6
ds <- DESeqDataSetFromMatrix(countData=counts,
  colData=covar, design=~ 0 + group)
cpm <- counts(ds) / rep(colSums(counts(ds)), each=nrow(counts(ds))) * 1e6
## take these genes where at least 3 samples have counts higher than 10
keep2 <- rowSums(cpm > 10) > 3
sum(keep2)
rs <- rowSums(counts(ds))
keep <- rs > 10
ds <- ds[ keep, ]
table(keep, keep2)
keep2 <- rowSums(cpm > 5) > 3
table(keep, keep2)
install.packages("pheatmap")
vst
foo <- rlog(counts)
foo <- rlog(ds)
head(foo)
counts(foo)
assayData(foo)
assays(foo)
?rlog
as.matrix(foo)
?DESeqTransform
assay(foo)
class(assay(foo))
head(assay(foo))
colSums(assay(foo))
foo2 <- assay(foo)
plot(foo2[,1], foo2[,2])
abline(0, 1)
resultsNames(ds)
ds <- DESeq(ds)
resultsNames(ds)
?result
?results
ares <- results(ds, contrast=c("group", "A549.RSV", "A549.mock"))
res <- results(ds, contrast=c("group", "A549.RSV", "A549.mock"))
head(res)
head(res[ order(res$pvalue), ])
res$gene_id <- rownames(res)
head(annotation)
res <- merge(annotation, res, by.x="PrimaryID", by.y="gene_id")
res <- merge(annotation, as.data.frame(res), by.x="PrimaryID", by.y="gene_id")
head(res)
head(res[ order(res$pvalue), ])
tmodCEres[ order(res$pvalue), ]
library(tmod)
tmodCERNOtest(res[ order(res$pvalue), ]$SYMBOL)
tmodCERNOtest(res[ sample(order(res$pvalue)), ]$SYMBOL)
?contrast
?results
savehistory()
zz()
nhbe <- counts[,1:6]
nhbe.cov <- covar[1:6,]
ds <- DESeqDataSetFromMatrix(nhbe, nhbe.cov, ~ group)
ds <- DESeq(ds)
resultsNames(ds)
res <- results(ds, name="group_NHBE.SC2V_vs_NHBE.mock")
head(res)
res$ID <- rownames(res)
res <- merge(annotation, res, by.x="PrimaryID", by.y="ID")
res <- merge(annotation, as.data.frame(res), by.x="PrimaryID", by.y="ID")
head(res)
res[ order(res$pvalue), ]
head(res[ order(res$pvalue), ])
tapply(nhbe["ENSG00000115009",], nhbe.cov$group, mean)
tapply(nhbe["ENSG00000115009",], nhbe.cov$group, mean, na.rm=TRUE)
nhbe[ "ENSG00000115009", ]
tapply(as.matrix(nhbe)["ENSG00000115009",], nhbe.cov$group, mean, na.rm=TRUE)
ds <- DESeqDataSetFromMatrix(nhbe, nhbe.cov, ~ 0 + group)
ds <- DESeq(ds)
resultsNames(ds)
res <- results(ds, name="groupNHBE.SC2V")
head(res)
res <- results(ds, contrast=c("group", "NHBE.SC2V", "NHBE.mock"))
res2 <- results(ds, contrast=c("group", "NHBE.mock", "NHBE.SC2V"))
cor(res$log2FoldChange, res2$log2FoldChange)
cor(res$log2FoldChange, res2$log2FoldChange, use="p")
?results
res2 <- results(ds, contrast=c(-1, 1))
resultsNames(ds)
cor(res$log2FoldChange, res2$log2FoldChange, use="p")
covar
ds <- DESeqDataSetFromMatrix(counts, covar, ~ group)
ds <- DESeq(ds)
resultsNames(ds)
resOnlyNHBE <- res
res1 <- results(ds, contrast=c(rep(0, 6), -1, 1))
cor(res1$log2FoldChange, resOnlyNHBE$log2FoldChange)
cor(res1$log2FoldChange, resOnlyNHBE$log2FoldChange, use="p")
ds <- DESeqDataSetFromMatrix(counts, covar, ~ 0 + group)
ds <- DESeq(ds)
resultsNames(ds)
res2 <- results(ds, c(rep(0, 6), -1, 1))
cor(res1$log2FoldChange, res2$log2FoldChange, use="p")
head(res1)
head(res2)
cor(res1$pvalue, res2$pvalue, method="s", use="p")
pal <- brewer.pal("Set2", n=8)
pal.t <- paste0(pal, '33')
norms <- lapply(1:3, function(x) rnorm(1000, sd=x))
norms.d <- lapply(norms, density)
plot(NULL, xlim=c(-10, 10), ylim=c(0, .5), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(norms.d[[x]], col=pal.t[x], border=NA)
  lines(norms.d[[x]], col=pal[x], lwd=3)
})
legend("topright", paste0("SD=", 1:3), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 200)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
alpha <- c(10, 100, 1000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 10000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 100000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 100000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 100000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(10, 100, 100000)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
alpha <- c(1, 10, 100)
mus <- rep(100, 3)
rnb2 <- lapply(1:3, function(x) rnbinom(1000, size=alpha[x], mu=mus[x]))
rnb.d2 <- lapply(rnb2, density)
plot(NULL, xlim=c(0, 300), ylim=c(0, .15), bty="n", 
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
legend("topright", sprintf("Mu=%d, alpha=%d", mus, alpha), col=pal, lwd=3, bty="n")
plot(NULL, xlim=c(0, 600), ylim=c(0, .03), bty="n",
  xlab="", ylab="density")
lapply(1:3, function(x) {
  polygon(rnb.d2[[x]], col=pal.t[x], border=NA)
  lines(rnb.d2[[x]], col=pal[x], lwd=3)
})
plotDispEsts(ds)
plotDispEsts(ds, bty="n")
##------ Tue May 12 21:13:23 2020 ------##
nhbe <- counts[, 1:6]
nhbe.cov <- covar[1:6, ]
nhbe <- counts[, 1:6]
nhbe.cov <- covar[1:6, ]
ds.nhbe <- DESeqDataSetFromMatrix(nhbe, nhbe.cov, ~ group)
resultsNames(ds.nhbe)
library(DESeq2)
nhbe <- counts[, 1:6]
nhbe.cov <- covar[1:6, ]
ds.nhbe <- DESeqDataSetFromMatrix(nhbe, nhbe.cov, ~ group)
resultsNames(ds.nhbe)
ds.nhbe <- DESeq(ds.nhbe)
resultsNames(ds.nhbe)
res1 <- results(ds.nhbe, name="group_NHBE.SC2V_vs_NHBE.mock")
xclip(1)
ds.nhbe <- DESeqDataSetFromMatrix(nhbe, nhbe.cov, ~ 0 + group)
ds.nhbe <- DESeq(ds.nhbe)
resultsNames(ds.nhbe)
res2 <- results(ds.nhbe, contrast=c("group", "NHBE.SC2V", "NHBE.mock"))
par(mfrow=c(1,2), bty="n")
plot(res1$log2FoldChange, res2$log2FoldChange, main="log2FoldChange")
plot(res1$pvalue, res2$pvalue, main="P-Value", log="xy")
abline(0, 1)
par(mfrow=c(1,2), bty="n", pch=19, col="#33333333")
plot(res1$log2FoldChange, res2$log2FoldChange, main="log2FoldChange")
abline(0, 1)
plot(res1$pvalue, res2$pvalue, main="P-Value", log="xy")
abline(0, 1)
?DESeq
?DESeqDataSetFromMatrix
design(ds)
design(ds) <- ~ group
design(ds) <- ~ group
ds <- DESeq(ds)
resultsNames(ds)
resultsNames(ds.nhbe)
resultsNames(ds.nhbe)
res3 <- results(ds.nhbe, contrast=c(-1, 1))
identical(res2, res3)
plot(res2$log2FoldChange, res3$log2FoldChange)
identical(as.data.frame(res2), as.data.frame(res3))
identical(as.data.frame(res1), as.data.frame(res3))
resultsNames(ds)
res4 <- results(ds, contrast=c(rep(6, 0), -1, 1))
length(resultsNames(ds))
c(rep(6, 0), -1, 1)
?results
## we can change the design on the fly
design(ds) <- ~ group
ds <- DESeq(ds)
resultsNames(ds)
res4 <- results(ds, contrast=c(rep(0, 6), -1, 1))
par(mfrow=c(1,2), bty="n", pch=19, col="#33333333")
plot(res1$log2FoldChange, res4$log2FoldChange, main="log2FoldChange")
abline(0, 1)
plot(res1$pvalue, res4$pvalue, main="P-Value", log="xy")
abline(0, 1)
par(mfrow=c(1,2), bty="n", pch=19, col="#33333333")
plot(res1$log2FoldChange, res4$log2FoldChange, main="log2FoldChange",
  xlab="Small data set", ylab="Large data set")
abline(0, 1)
plot(res1$pvalue, res4$pvalue, main="P-Value", log="xy",
  xlab="Small data set", ylab="Large data set")
abline(0, 1)
x <- factor(rep(c("A", "B"), each=5))
y <- rep(c(1,3), each=5) + rnorm(10)
anova(lm(y ~ x))
summary(lm(y ~ x))
aov(y ~ x)
summary(aov(y ~ x))
anova(lm(y ~ 0 + x))
summary(lm(y ~ x))
t.test(y[1:5], y[6:10])
summary(lm(y ~ x))
mean(y[1:5])
mean(y[6:10])
mean(y[6:10])
mean(y[6:10])-mean(y[1:5])
x <- factor(rep(c("A", "B"), each=5))
y <- rep(c(4,6), each=5) + rnorm(10)
anova(lm(y ~ x)) 
set.seed(123)
x <- factor(rep(c("A", "B"), each=5))
y <- rep(c(4,6), each=5) + rnorm(10)
anova(lm(y ~ x)) 
set.seed(123)
x <- factor(rep(c("A", "B"), each=5))
y <- rep(c(4, 7), each=5) + rnorm(10)
anova(lm(y ~ x))
t.test(y[1:5], y[6:10])
summary(lm(y ~ x))
summary(lm(y ~ 0 + x))
contrast(lm(y ~ 0 + x), c(-1, 1))
contrasts(lm(y ~ 0 + x), c(-1, 1))
summary(lm(y ~ 0 + x))
anova(lm(y ~ 0 + x))
res <- results(ds, contrast=c("group", "A549.RSV", "A549.mock"))
rmarkdown::render("../Radvanced/Lectures/Lecture_02/lecture_02.rmd")
res4
head(res4[ order(res4$pvalue), ])
ds <- DESeqDataSetFromMatrix(countData=counts, 
  colData=covar, design=~ 0 + group)
ds <- DESeq(ds)
resultsNames(ds)
foo <- results(ds, contrast=c("group", "NHBE.SC2V", "NHBE.mock"))
head(foo[ order(foo$pvalue), ])
