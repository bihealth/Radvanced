##------ Wed May 20 13:11:44 2020 ------##

x <- rnorm(10)
sqrt(sum((x - mean(x))^2)/length(x))
sd(x)
sqrt(sum((x - mean(x))^2)/(length(x)-1))
options(digits=4)
sqrt(sum((x - mean(x))^2)/(length(x)-1))
sd(x)
sqrt(sum((x - mean(x))^2)/length(x))
    sd_pop <- function(x, ...) {
      n <- length(x)
      return(sqrt((n-1)/n)*sd(x))
    }
sd_pop(x)
timestamp()
savehistory()
foo <- replicate(100, { rnorm(10, sd=1, mean=0) }, simplify)
foo <- replicate(100, { rnorm(10, sd=1, mean=0) }, simplify=FALSE)
class(foo)
length(foo)
foo[[1]]
sd_pop(foo[[1]])
sd(foo[[1]])
samps <- foo
      sd.p <- sapply(samps, sd_pop)
      sd.s <- sapply(samps, sd)
head(sd.p)
head(sd.s)
    mysamp <- function(N, rep=1000, mean=0, sd=1) {
      ## you could do the following with lapply, but replicate is more
      ## elegant
      samps <- replicate(rep, { rnorm(N, mean=mean, sd=sd) }, simplify=FALSE)
      sd.p <- sapply(samps, sd_pop)
      sd.s <- sapply(samps, sd)
      return(data.frame(SD_pop=mean(sd.p), SD_samp=mean(sd.s)))
    }
blah <- mysamp(10, rep=3)
blah
    mysamp <- function(N, rep=1000, mean=0, sd=1) {
      ## you could do the following with lapply, but replicate is more
      ## elegant
      samps <- replicate(rep, { rnorm(N, mean=mean, sd=sd) }, simplify=FALSE)
      sd.p <- sapply(samps, sd_pop)
      sd.s <- sapply(samps, sd)
      return(data.frame(SD_pop=mean(sd.p), SD_samp=mean(sd.s)))
    }
    N <- 3:200
    res <- sapply(N, mysamp)
dim(res)
head(res)
    plot(N, res["SD_samp", ])
    points(N, res["SD_pop", ], col="red")
    abline(h=1)
    legend("bottomright", c("SD_samp", "SD_pop"), col=c("black", "red"), pch=1)
x <- rnorm(100)
quantile(x, .75)
quantile(x, .25)
x <- quantile(x, .25)
x <- rnorm(100)
sum(x < quantile(x, .25))
quantile(x, c(.25, .75))
x <- log(rnorm(100))
x <- log(rnorm(100) + 100)
plot(density(x))
x <- log(rnorm(100) + 10)
plot(density(x))
x <- log(rnorm(100) + 5)
plot(density(x))
x <- log(rnorm(100, mean=100, sd=20) + 5)
plot(density(x))
median(x)
mean(x)
mean(rnorm(100, mean=100, sd=20))
median(rnorm(100, mean=100, sd=20))
x <- 2^(rnorm(100, mean=100, sd=20) + 5)
plot(density(x))
x <- 2^(rnorm(100) + 5)
plot(density(x))
abline(v=mean(x))
abline(v=median(x))
ds
res <- results(ds, contrast=c("group", "A549.RSV", "A549.mock"))
res <- merge(annot, data.frame(res), by.x="PrimaryID", by.y=0)
resultsNames(res)
resultsNames(ds)
matrix(0, nrow=8, ncol=4)
cmtx <- matrix(0, nrow=8, ncol=4)
rownames(cmtx) <- resultsNames(ds)
colnames(cmtx) <- c("NHBE.SC2V", "A549.SC2V", "A549.RSV", "A549.IAV")
cmtx
cmtx[,1] <- c(rep(0, 6), -1, 1)
cmtx
cmtx[,2] <- c(0, 0, -1, 0, 0, 1, 0, 0)
cmtx[,3] <- c(0, 0, 0, 1, -1, 0, 0, 0)
cmtx[,4] <- c(1, -1, rep(0, 6))
cmtx
res_full <- apply(cmtx, 2, function(contrast) results(ds, contrast=contrast))
class(res_full)
names(res_full)
res_full$NHBE.SC2V
results(ds, contrast=list("groupA549.IAV", "groupA549.IAV.mock"))
resultsNames(ds)
cmtx
results(ds, contrast=c(1, -1, 0, 0, 0, 0, 0, 0))
log(100/10)
log(1000/100)
res_full
names(foo)
names(res_full)
all(rownames(res_full[[1]]) == rownames(res_full[[2]]))
foo <- data.frame(NHBE.lfc=res_full$NHBE.SC2V$log2FoldChange, NHBE.pval=res_full$NHBE.SC2V$padj, A549.lfc=res_full$A549.SC2V$log2FoldChange, res_full$A549.SC2V$padj)
head(foo)
rownames(foo) <- rownames(res_full[[1]])
head(foo)
foo.CCC <- foo[ abs(foo$NHBE.lfc) < .2 & abs(foo$A549.lfc) < .2, ]
head(foo.CCC)
lcounts <- rlog(counts(ds))
head(lcounts)
rownames(lcounts) <- rownames(ds)
head(lcounts)
showgene(lcounts["ENSG00000000460", ], group=covar$group)
par(mar=c(7,3,3,3)
par(mar=c(7,3,3,3))
showgene(lcounts["ENSG00000000460", ], group=covar$group)
foo.AAA <- foo[ abs(foo$NHBE.lfc - 2) < .2 & abs(foo$A549.lfc - 2 < .2), ]
dim(foo.AAA)
head(foo.AAA)
tail(foo.AAA)
foo.AAA <- foo[ abs(foo$NHBE.lfc - 2) < .5 & abs(foo$A549.lfc - 2 < .5), ]
foo.AAA
foo.AAA <- foo.AAA[ !is.na(foo.AAA[,1]), ]
foo.AAA
head(foo.AAA[ foo.AAA$NHBE.pval < .01, ])
foo.AAA <- foo[ abs(foo$NHBE.lfc - 2) < .5 & abs(foo$A549.lfc - 2) < .5, ]
foo.AAA <- foo.AAA[ !is.na(foo.AAA[,1]), ]
head(foo.AAA)

foo.DDD <- foo[ abs(foo$NHBE.lfc) < .2 & abs(foo$A549.lfc - 2) < .5, ]
foo.DDD <- foo.DDD[ !is.na(foo.DDD[,1]), ]
foo.DDD
head(foo.DDD)
foo <- data.frame(NHBE.lfc=res_full$NHBE.SC2V$log2FoldChange, NHBE.pval=res_full$NHBE.SC2V$padj, A549.lfc=res_full$A549.SC2V$log2FoldChange, A549.pval=res_full$A549.SC2V$padj)
foo.DDD <- foo[ abs(foo$NHBE.lfc) < .2 & abs(foo$A549.lfc - 2) < .5, ]
foo.DDD <- foo.DDD[ !is.na(foo.DDD[,1]), ]
head(foo.DDD[ order(foo.DDD$A549.pval), ])
rownames(foo) <- rownames(res_full[[1]])
foo.DDD <- foo[ abs(foo$NHBE.lfc) < .2 & abs(foo$A549.lfc - 2) < .5, ]
foo.DDD <- foo.DDD[ !is.na(foo.DDD[,1]), ]
head(foo.DDD[ order(foo.DDD$A549.pval), ])
covar
group <- covar$group
group <- factor(group, levels=unique(sort(group)))
group
group <- factor(group, levels=unique(group))
group
head(foo.AAA)
head(foo.DDD[ order(foo.DDD$A549.pval), ])
showgene(lcounts[ "ENSG00000137628", ], group=group)
showgene(lcounts[ "ENSG00000149131", ], group=group)
names(res_full)
resultsNames(ds)
groupNHBE.SC2V-groupNHBE.mock - groupA549.SC2V +groupA549.mock
resultsNames(ds)
resultsNames(ds)
contrast <- c(0, 0, 1, 0, 0, -1, -1, 1)
res.interaction <- results(ds, contrast=contrast)
head(res.interaction[ order(res.interaction$pvalue), ])
showgene(lcounts[ "ENSG00000187608", ], group=group)
showgene(lcounts[ "ENSG00000115009", ], group=group)
showgene(lcounts[ "ENSG00000163735", ], group=group)
head(foo)
interesting <- sign(foo$NHBE.lfc) != sign(foo$A549.lfc)
head(foo[ interesting, ])
interesting <- (sign(foo$NHBE.lfc) != sign(foo$A549.lfc)) & !is.na(foo$NHBE.lfc)
head(foo[ interesting, ])
all(rownames(foo) == rownames(res.interaction))
head(res.interaction[ res.interaction$padj < .05 & interesting, ])
which(is.na(interesting))
head(res.interaction[ !is.na(res.interaction$padj) & res.interaction$padj < .05 & interesting, ])
discordant <- res.interaction[ !is.na(res.interaction$padj) & res.interaction$padj < .05 & interesting, ]
discordant[ order(discordant$padj), ]
showgene(lcounts[ "ENSG00000134070", ], group=group)
showgene(lcounts[ "ENSG00000145901", ], group=group)
showgene(lcounts[ "ENSG00000170477", ], group=group)
annot[ "ENSG00000170477", ]
annot[ annot$PrimaryID == "ENSG00000170477", ]
showgene(lcounts[ "ENSG00000128342", ], group=group)
discordant <- merge(annot, discordant, by.x="PrimaryID", by.y=0)
discordant <- merge(annot, as.data.frame(discordant), by.x="PrimaryID", by.y=0)
head(discordant)
tmodHGtest(discordant$SYMBOL, annot$SYMBOL)
discordant <- res.interaction[ !is.na(res.interaction$padj) & res.interaction$padj < .15 & interesting, ]
discordant <- merge(annot, as.data.frame(discordant), by.x="PrimaryID", by.y=0)
tmodHGtest(discordant$SYMBOL, annot$SYMBOL)
tmodHGtest(fg=discordant$SYMBOL, bg=annot$SYMBOL)
tmodHGtest(fg=discordant$SYMBOL, bg=setdiff(annot$SYMBOL, discordant$SYMBOL))
