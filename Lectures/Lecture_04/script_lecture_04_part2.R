##------ Fri May 22 10:00:45 2020 ------##

## --------------------------------------------------------------------------
## If you are missing any objects or commands, please run the code from the
## lecture! Also the parts which are hidden, meaning that you need to download
## the file Lecture_04.rmd, open it in your rstudio session  and rund the 
## code chunks.
## --------------------------------------------------------------------------

library(tmod)
cgs
tcgs
tcgs$MODULES
tcgs$MODULES2GENES
tcgs$GENES
tcgs$GENES2MODULES
 head(l)
tcgs$MODULES2GENES$Jak_Stat
l <- res$SYMBOL[ order(res$padj) ]
library(tmod)
tmodCERNOtest(l, mset=tcgs)
library(msigdbr)
h.df <- msigdbr("Homo sapiens", "H")
dim(h.df)
head(h.df)
length(unique(h.df$gs_name))
head(h.cgs)
h.cgs <- lapply(h.names, function(n) 
h.df$human_gene_symbol[ h.df$gs_name == n ])
h.cgs <- lapply(h.names, function(n) 
head(h.cgs)
h.names <- unique(h.df$gs_name)
h.cgs <- lapply(h.names, function(n) 
  h.df$human_gene_symbol[ h.df$gs_name == n ])
head(h.cgs)
head(h.names)
h.names <- gsub("HALLMARK_", "", h.names)
names(h.cgs) <- h.names
htmod <- makeTmod(modules=data.frame(ID=h.names, Title=h.names), modules2genes=h.cgs)
htmod
dim(h.df)
htmod
htmod[ 1:10 ]
head(cgs.df)
head(res.cphg)
vec <- nrow(res):1
names(vec) <- res$SYMBOL[ order(res$pvalue) ]
res.cp <- GSEA(vec, TERM2GENE=cgs.df)
library(clusterProfiler)
vec <- nrow(res):1
names(vec) <- res$SYMBOL[ order(res$pvalue) ]
res.cp <- GSEA(vec, TERM2GENE=cgs.df)
cgs
res.cp
head(res.cp)
head(res.cp)[,1:5]
head(res.cp)[,1:7]
vec <- -log10(res$pvalue)
names(vec) <- res$SYMBOL
vec <- sort(-vec)
head(vec)
vec <- sort(vec, decreasing=T)
head(vec)
vec <- -log10(res$pvalue)
head(vec)
names(vec) <- res$SYMBOL
vec[is.na(vec)] <- -log10(1)
head(vec)
vec <- sort(vec, decreasing=T)
head(vec)
vec[is.infinite(vec)] <- 290.9112
head(vec)
res.cp <- GSEA(vec, TERM2GENE=cgs.df)
head(res.cp)
head(all_eig)
group1 <- covar$group == "RSV.mock"
sum(group1)
covar$group
group1 <- covar$group == "A549.RSV.mock"
group2 <- covar$group == "A549.RSV"
foo <- apply(all_eig, 1, function(x) t.test(x[group1], x[group2])$p.value)
foo
?tmodPLAGEtest
res_tmod <- tmodCERNOtest(l)
res_tmod
m1 <- "LI.M4.1"  # cell cycle (I)
m2 <- "LI.M4.2"  # PLK1 signaling events
m3 <- "LI.M22.0" # mismatch repair
glist <- getModuleMembers(c(m1, m2, m3))
intersect(glist[[1]], glist[[2]])
res[ res$SYMBOL %in% intersect(glist[[1]], glist[[2]]), ]
library(VennDiagram)
## grid graphics is between basic R and ggplot.
## No need to understand it.
grid.newpage()
venn.plot <- venn.diagram(glist, NULL,
  fill=brewer.pal(3, "Pastel1"),
  cex=1.5, cat.cex=1.5, 
  lty="blank",
  fontfamily="sans", cat.fontfamily="sans")
library(RColorBrewer)
library(VennDiagram)
## grid graphics is between basic R and ggplot.
## No need to understand it.
grid.newpage()
venn.plot <- venn.diagram(glist, NULL,
  fill=brewer.pal(3, "Pastel1"),
  cex=1.5, cat.cex=1.5, 
  lty="blank",
  fontfamily="sans", cat.fontfamily="sans")
grid.draw(venn.plot)
simmat
simmat[1:10,1:10]
simmat[1:10,1:10]
groups
tmod_res[ tmod_res$ID %in% groups$LI.M72.0, ]
res_tmod[ res_tmod$ID %in% groups$LI.M72.0, ]
res_tmod[ res_tmod$ID %in% groups$LI.M43.1, ]
res_tmod[ res_tmod$ID %in% groups$LI.M4.0, ]
upset(res_tmod$ID[1:15], min.overlap=5, group=TRUE)
upset(res_tmod$ID[1:15], min.overlap=5, group=TRUE, value="j")
upset(res_tmod$ID[1:15], min.overlap=5, group=TRUE, value="j", cutoff=.1)
upset(res_tmod$ID, min.overlap=5, group=TRUE, value="j", cutoff=.1)
upset(res_tmod$ID, min.overlap=5, group=TRUE, value="j", cutoff=.1, max.comb=3)
upset(res_tmod$ID, min.overlap=5, group=TRUE, value="j", cutoff=.1, max.comb=3, min.group=2)
rmarkdown::render("lecture_04.rmd")
rmarkdown::render("lecture_04.rmd")
rmarkdown::render("lecture_04.rmd")
rmarkdown::render("lecture_04.rmd")
upset(lea, min.overlap=5, group=TRUE, labels=res_tmod$Title[1:25], value="j")
savehistory()
