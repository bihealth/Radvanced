Before you start, run the following commands to make sure we are on the
same page. Essentially, this is the DESeq2 analysis of the Covid data.

```
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
 })

## get the newest tmod
library(devtools)
install_github("https://github.com/january3/tmod")
```

 1. OK, first a couple of questions regarding the code above.

     1. What is the `res` object`? How many elements does it have? 
     2. Can you say which five genes have the lowest p-value when
        comparing IAV with mock treatment in A549 cells?
     3. Which genes are significant in both NHBE and A549 cells when
        treated with SC2V?
     3. Plot the log2FoldChange in SC2 treatment of A549 vs NHBE cells.
        Mark the genes which are significant in *both* comparisons with
        red.

 2. Basic gene set enrichment with tmod.

     1. First, run the tmodCEROtest on the list of genes ordered by p-value
        in the SC2V treatment of A549 cells. Which groups of genes are
        significant?

     2. Choose one of the gene sets to study further (better choose one
        with not too many genes). Which genes belong to the gene set? (use
        `getModuleMembers`). Which of these genes have a significant
        p-value?

     3. Plot the ROC curve (function `evidencePlot`). The first parameter
        is the ordered gene list, the second, the gene set (module) ID.
        With the parameter `gene.labels=TRUE` you can show the gene
        position on the ordered list.

     4. You can color the genes using the parameter `gene.color`. This
        parameter takes a named vector, where colors are the values and
        genes are the names. Create a vector such that genes which are
        significant are colored in red, and genes which are not are grey
        and make an evidence plot.

 3. Finding groups of gene sets

     1. Use the function `modGroups` the gene sets obtained in the previous
        step. Try different values of the `min.overlap` parameter. What
        groups of gene sets do you see? Can you show what the gene sets do?
        (Hint: if `tmod.res` is the data frame with tmod results, then 
        `tmod.res[ c("id1", "id2"), ]$Title` shows the description of the
        gene sets with the IDs id1 and id2.

     2. Use the `upset` function to generate an upset plot. This function
        has many parameters. Try changing the `cutoff` parameter and the
        `min.overlap` parameter. Can you interpret what is on the plot?

 4. tmod Panel plots (function `tmodPanelPlot`) are an important tool to
    compare the enrichments in different contrasts. This part of the
    homework will show you how to generate them.

     1. Using `lapply`, generate the enrichment for every contrast in the
        `res` list. The result should be a list called `tmod_res`, and each
        element of the list should be the results of tmod gene set
        enrichment. This will be the main object passed to `tmodPanelPlot`.

     2. The basic usage of `tmodPanelPlot` is `tmodPanelPlot(tmod_res)`,
        where `tmod_res` is a list of tmod results. However, there are
        several parameters to tweak the output. Please play with the
        following:

          1. Adapt the font size with `text.cex` (best use values below
             1.0, e.g. 0.5)
          2. What does `grid="b"` do?
          3. `pval.thr` determines the maximum p-value which will be shown
             on the plot. Change it to more conservative (lower) values.
          4. `filter.rows.pval` and `filter.rows.auc` are used to remove
             gene sets from the plot which do not fullfill the specified
             criteria, i.e. p value lower than the threshold in one or more
             of the comparisons shown, and AUC above the threshold in one
             or more of the comparisons. Modify the values to show only the
             "best" enrichments.

      3. To show whether genes in a gene set go up or down, we need to
         provide this information to tmod. We do this contructing an object
         called `pie` which holds the information about how many genes in a
         gene significantly set go up or down. To construct the object, do
         the following:

           1. Create a matrix calleds `lfcs` with log fold changes from
              each contrast. Use the `sapply` function, for example 
              `sapply(res, function(x) x$log2FoldChange)`.

           2. Similarly, construct an object containing adjusted p-values
              (call it `pvals`). Make sure both lfcs and pvals do not
              contain any NA's! For log2FoldChange, replace the NA by 0;
              for p-vals, replace the NA's with 1. (for example, do 
              pvals[is.na(pvals)] <- 1`).

           3. Now call `pie <- tmodDecideTests(res[[1]]$SYMBOL, lfcs, pvals)`

           4. Take a look at the resulting object. It is a list. What are
              the members of that list? Can you test the results manually?

           5. Call `tmodPanelPlot(tmod_res, pie=pie, grid="b")`. Try to
              tune the parameters to get a nice plotting result.

 4. I am aware that this is a lot of work, but trust me, it is worth it:

      1. Add another contrast to `res`: the interaction between SC2V
         infection and cell type.

      1. Create a tmod database from MSigDB which contains only the
         REACTOME gene sets.
      
      2. Repeat *all* the previous steps for this new set of gene sets and
         contrasts using the `mset` parameter. Are the results similar?
         What would be your conclusion about how cells react to the
         infection?


    





