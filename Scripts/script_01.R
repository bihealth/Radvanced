##------ Wed May  6 10:18:58 2020 ------##
## read the counts and annotation
counts <- read.table("Data/covid_data/counts.tsv")
annotation <- read.csv("Data/covid_data/annotation.all.csv")

## check that the two data frames are in line
all(rownames(counts) == annotation$PrimaryID)

## otherwise we use "match" to align them
head(annotation$SYMBOL[ match(rownames(counts), annotation$PrimaryID)])

## OK, so you did not have the covariate table. No worries, we are making
## it up on the spot.
cell.types <- c(rep("NHBE", 6), rep("A549", 14))
treatment <- c(rep(c("M", "S2"), each=3), rep(c("M", "S2"), each=3), rep(c("M", "RSV"), each=2), rep(c("M", "IAV"), each=2))
covariate_table <- data.frame(ID=colnames(counts), cell.type=cell.types, treatment=treatment)

## you can also take a look at the original covariate table. 
## file Data/covid_data/covariate_file.txt

## let us take a look at the density plots.
plot(density(counts[,1]))
plot(density(log2(counts[,1])))
lines(density(log2(counts[,2])), col="red")
lines(density(log2(counts[,3])), col="blue")
lines(density(log2(counts[,4])), col="green")

## we could automatise it
d <- apply(log2(counts), 2, density)
plot(NULL, xlim=c(0, max(log2(counts))), ylim=c(0, 0.2))
lapply(1:20, function(i) lines(d[[i]], col=i))

## ----- ggplot version
## for ggplot, we need a data frame. Preparing it is sometimes pretty
## awful.

## Basically, we have two ways. One is calculating density with density()
## and then drawing the lines. We will use the d list for that.

d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=d[[i]]$x, y=d[[i]]$y, ID=colnames(counts)[i]))
## please take a look at the resulting objects, i.e. for example d.gg[[1]]

## fuse it into one data frame
d.gg <- Reduce(rbind, d.gg)

## now we can plot
ggplot(d.gg, aes(x=x, y=y, group=ID)) + geom_line(aes(col=ID))

## alternatively, we convert the counts matrix into a data frame
## and use geom_density
d.gg <- lapply(1:ncol(counts), function(i) data.frame(x=counts[,i], ID=colnames(counts)[i]))
d.gg <- Reduce(rbind, d.gg)
d.gg$x <- log2(d.gg$x)
## check how large this data frame is!
dim(d.gg) 

ggplot(d.gg, aes(x=x, group=ID)) + geom_density(aes(col=ID))

## ----- Let us look at a contrast
foo <- read.table("Data/covid_data/contrast.A549_RSV_ID2.csv", header=TRUE)

## add gene description to the contrast results
foo$Description <- annotation$GENENAME[ match(foo$gene_id, annotation$PrimaryID)]

## we could also use merge. Note that the matching columns have different
## names!
foo.m <- merge(foo, annotation, by.x="gene_id", by.y="PrimaryID")

## counts has been read by read.table, so it is a data frame. For many
## purposes it is easier to use it as a matrix, though
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

## how about plotting by cell type *and* treatment?
covariate_table$c.t <- paste0(covariate_table$cell.type, ".", covariate_table$treatment)

## but actually, we should log the data!
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y")

## we can add a beeswarm
library(beeswarm)
beeswarm(unlist(counts[id, ]) ~ covariate_table$c.t, log=T, add=T)

## make a nicer version

boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, ylab="Raw counts", lab.cex=.7)
## make a nicer version
library(RColorBrewer)
boxplot(mtx[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, xlab="Experimental group", ylab="Raw counts", cex.axis=.7)
beeswarm(mtx[id, ] ~ covariate_table$c.t, log=T, add=T, 
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)

## ggplot version
library(ggbeeswarm) # install it if you don't have it
library(cowplot) # to clean up the output
gg <- data.frame(covariate_table, y=mtx[id,])
ggplot(gg, aes(x=c.t, y=y)) + geom_boxplot(aes(col=c.t)) +
  geom_beeswarm(aes(col=c.t), size=2.5) +
  scale_y_continuous(trans=log2_trans())

## actually, these are *incorrect* plots!
## they show the RAW counts. However, count strongly depend on the library
## size (total number of reads in a sample).
counts2 <- read.table("Data/covid_data/DESeq2.all.rld.model.csv", header=TRUE, row.names=1)
mtx2 <- data.matrix(counts2[, 1:20])

boxplot(mtx2[id, ] ~ covariate_table$c.t, log="y", frame.plot=FALSE, 
  main=id, xlab="Experimental group", ylab="Raw counts", cex.axis=.7)
beeswarm(mtx2[id, ] ~ covariate_table$c.t, log=T, add=T, 
  col=brewer.pal(6, "Dark2"), pch=19, cex=1.4)


