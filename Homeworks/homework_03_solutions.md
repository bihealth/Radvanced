 1. **Answer:** The `sd()` function returns the sample SD defined as 

    <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\frac{\sum_i=1^n(x_i-\overline{x})^2}{n - 1}}">

    To easily calculate the population SD using `sd()`, note that

    <img src="https://render.githubusercontent.com/render/math?math=\text{SD}_\text{population}=\sqrt{\frac{\sum_{i=1}^n(x_i-\overline{x})^2}{n}}=\sqrt{\frac{n-1}{n-1}\cdot\frac{\sum_{i=1}^n(x_i-\overline{x})^2}{n}}=\sqrt{\frac{n-1}{n}}\cdot\sqrt{\frac{\sum_{i=1}^n(x_i-\overline{x})^2}{n-1}}=\sqrt{\frac{n-1}{n}}\cdot\text{SD}_\text{sample}">

    Therefore, to calculate population SD, simply use
    `sqrt((n-1)/n)*sd(x)`. We can code it as a function:

    ```
    sd_pop <- function(x, ...) {
      n <- length(x)
      return(sqrt((n-1)/n)*sd(x))
    }
    ```

 1. **Answer:** Here is one possible algorithm:
    
     1. For a given N, sample the normal distribution a hundred times (with known SD as a
        parameter; use the `rnorm` function), each time drawing N numbers.
        
     2. Calculate the SD with both methods for each of the hundred samples
        and calculate the average for both estimators.

     2. Plot the averages of both estimators against N. Draw a horizontal
        line (e.g. with `abline`) at the SD that was used to generate the
        samples.

    ```
    ## define a function for sampling
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
    dim(res) ## two rows, 198 columns!
    par(bty="n")
    plot(N, res["SD_samp", ])
    points(N, res["SD_pop", ], col="red")
    abline(h=1)
    legend("bottomright", c("SD_samp", "SD_pop"), col=c("black", "red"), pch=1)
    ```

    As you can see, the population based estimate (with "N" instead of
    "N-1") consistently underestimates the real variance.
    
 2. Using the approach presented in the lecture and the *full* data set,
    create the object `ds` (DESeq object) with design `~ 0 + group`.
    Create results for the following contrasts and store the results in a
    list:

     * NHBE.SC2V vs NHBE.mock
     * A549.SC2V vs A549.mock
     * A549.RSV vs A549.RSV.mock
     * A549.IAV vs A549.IAV.mock

    **Answer:** There are many ways to do that. I like the approach with a
    contrast matrix as below:

    ```
    coefs <- resultsNames(ds)
    cmtx <- matrix(0, nrow=length(coefs), ncol=4)
    rownames(cmtx) <- resultsNames(ds)
    colnames(cmtx) <- c("NHBE.SC2V", "A549.SC2V", "A549.RSV", "A549.IAV")
    cmtx[,1] <- c(rep(0, 6), -1, 1)        # +groupNHBE.SC2V -groupNHBE.mock
    cmtx[,2] <- c(0, 0, -1, 0, 0, 1, 0, 0) # +groupA549.SC2V -groupA549.mock
    cmtx[,3] <- c(0, 0, 0, 1, -1, 0, 0, 0) # +groupA549.RSV  -groupA549.RSV.mock
    cmtx[,4] <- c(1, -1, rep(0, 6))        # +groupA549.IAV  -groupA549.IAV.mock

    res_full <- apply(cmtx, 2, function(contrast) results(ds, contrast=contrast))
    ```
     
    **More info than you need:** You could use the function `makeContrasts`
    from the `limma` package to do that semi-automatically. Here is how:

    ```
    design <- model.matrix(~ 0 + group, data=covar)  ## check the output!
    all(colnames(design) %in% resultsNames(ds)) ## must be TRUE
    cmtx2 <- limma::makeContrasts(
                          NHBE.SC2V="groupNHBE.SC2V-groupNHBE.mock",
                          A549.SC2V="groupA549.SC2V-groupA549.mock",
                          A549.RSV="groupA549.RSV-groupA549.RSV.mock",
                          A549.IAV="groupA549.IAV-groupA549.IAV.mock",
                          levels=design)
    ## the cmtx2 is identical to cmtx
    ```


 2. **Question:** which of the genes show an interaction?

    **Answer:** The last three rows.

 3. Simply by filtering the results of direct comparisons between SC2V and
    mock for A549 and NHBE cells, try to find genes that roughly follow
    the six patterns in the table above. It does not need to be
    significant, and it does not have to be very precise. Create dot plots
    (like with `beeswarm`) for the selected genes.

    **Answer:**

    Here is one possible way of doing that. We create a table with the log2
    fold changes and p-values from the two contrasts corresponding to the
    given comparisons.

    ```
    cmptab <- data.frame(
      NHBE.lfc=res_full$NHBE.SC2V$log2FoldChange, 
      NHBE.pval=res_full$NHBE.SC2V$padj, 
      A549.lfc=res_full$A549.SC2V$log2FoldChange, 
      A549.pval=res_full$A549.SC2V$padj)
    rownames(cmptab) <- rownames(ds)
    ## remove all NAs
    cmptab <- cmptab[ complete.cases(cmptab), ]
    ```

    We can now create the different subsets. We can use `subset` to
    simplify the syntax somehow and `showgene` from `tmod` package for
    plotting:

    ```
    ## CCC: similar change in both contrasts
    ## search for genes which have in both contrasts lfc ~ 0
    subset.CCC <- subset(cmptab, abs(NHBE.lfc) < .2 & abs(A549.lfc) < .2)
    
    ## search for genes which have in both contrasts lfc ~ 2
    subset.AAA <- subset(cmptab, abs(NHBE.lfc - 2) < .5 & abs(A549.lfc - 2) < .5)

    ## genes with almost no change in NHBE but a large change in A549
    subset.DDD <- subset(cmptab, abs(NHBE.lfc) < .2 & abs(A549.lfc) > 1.5)
    head(subset.DDD[ order(subset.DDD$A549.pval), ])

    ## etc.
    ```

    For display, we need the normalized counts:

    ```
    lcounts <- rlog(counts(ds))
    rownames(lcounts) <- rownames(ds)

    group <- factor(covar$group, levels=unique(covar$group))

    showgene(lcounts["ENSG00000089127", ], group)
    ```
    



 4. Inspect the obtained results for an interaction test. Plot some of the
    top genes. Do you see the interaction? Look up the gene examples you
    have selected before in the results tables. What do the associated
    log2FoldChanges mean?


    **Answer:** See the script file from the lecture.
