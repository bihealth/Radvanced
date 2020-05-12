
 1. Load all the contrast files. Can you do it automatically? Look up the
    function `list.files()` and its parameter `pattern`. Use `lapply` to
    load all the files.

    **Solution:**

    ```
    files.list <- list.files("Data/covid_data/", pattern="contrast*", full.names=TRUE)
    files <- lapply(files.list, read.table, header=TRUE)
    ```

    **Explanations:**

    We need the `full.names` parameter to `list.files`, because otherwise
    the file names do not contain the path `Data/covid_data/`.

    Remember that `lapply` applies a function to each element of the list.
    Any additional arguments (here `header=TRUE`) are passed to that
    function.


 2. Create a matrix containing only the log2 fold changes from all
    contrasts.  Either use `lapply` and `Reduce`, or try `sapply` which can
    do some automatic transformations of the output.

    **Solution:**

    ```
    # reorder the tables
    genes <- files[[1]]$gene_id
    files <- lapply(files, function(x) x[ match(genes, x$gene_id), ])

    # extract log2FoldChange
    data <- lapply(files, function(x) x$log2FoldChange)
    # or
    data <- lapply(files, function(x) x[["log2FoldChange"]])
    # or
    data <- lapply(files, function(x) x[ , "log2FoldChange"]])
    # or
    data <- lapply(files, `[[`, "log2FoldChange")

    # join the tables
    mtx <- Reduce(cbind, data)

    ## alternatively
    mtx <- sapply(files, `[[`, "log2FoldChange")

    ## did it work?
    class(mtx)
    dim(mtx)
    ```

    **Explanations:**

    First thing to do is to order all the data frames in the same way. The
    problem is that the contrast results are sorted by log2FoldChange, but
    we need to have in the resulting matrix each row corresponding to one
    gene. Hence the first two statements. We could use `merge`, but it
    would have been more messy.

    Remember that a data frame (as returned by `read.table`) is a special
    kind of list, so you can access columns using `$`, `[[` and `[, ]`.
    Hence the first three expressions are equivalent.

    Remeber also that `[[` is a function, albeit using special characters
    as name. Therefore, you need to use back ticks to call it as a
    function, for example

    ```
    ## return the first element of files
    `[[`(files, 1)
    ```

    Thus, we can also simply pass it to `lapply` (but using backticks!) and
    give the "log2FoldChange" as parameter.

    Finally, `Reduce` takes as parameters a function and a list. Given a
    list with elements a, b and c, it applies the function to a and b and
    generates a new result, let us call it ab. Then, it applies the
    function to ab and c. And so forth. So in our case, it takes the first two
    vectors, applies `cbind` to them, and then applies cbind to the
    resulting matrix and the third vector, producing a three-column matrix.

    `sapply` is a convenient function, but since it converts the data
    automatically, we are not guaranteed that the result will be a matrix.
    Thus, we need to be wary and check the resulting type.

 3. Plot the correlation between two contrasts. Can you make a plot that
    shows all the correlations in one go? Either take a look at the web
    site
    [https://www.r-graph-gallery.com/correlogram.html](https://www.r-graph-gallery.com/correlogram.html)
    or (which is harder, but more flexible) use the `pairs` or `ggpairs`
    functions (`ggpairs` is in the GGally package).

    **Solution:**

    ```
    plot(mtx[,1], mtx[,2], pch=19)
    plot(mtx[,1], mtx[,4], pch=19)
    smoothScatter(mtx[,1], mtx[,4])

    ## not part of the homework, but important!
    cor(mtx[,1], mtx[,2], use="p")
    cor(mtx, use="p")
    cor(mtx, use="p", methods="spearman")

    ## simplest correlogram ever!
    ## #33333333 is transparent grey
    plot(as.data.frame(mtx), pch=19, col="#33333333")

    ## correlograms with packages
    library(corrgram)
    corrgram(mtx, lower.panel=panel.cor, upper.panel=panel.ellipse, diag.panel=panel.density)

    ## GGally
    library(GGally)
    mtx.df <- data.frame(mtx)
    ggpairs(mtx)

    ## using pairs
    pairs(mtx)

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

    ```

    **Explanations:**

    The `cor` function can take two arguments (vectors), or even a whole
    matrix. We need to use the `use="p"` parameter, which stands for
    "pairwise.complete.obs", because some of the rows are full of NA's.

    Correlogram are a great way to get an overview of the data. The figure
    is divided in panels. On the diagonal, there are the variable names
    (names of columns in the matrix, if existing). Each panel results from
    the combination (for example, a scatter plot) of two variables on the 
    diagonal. 

    There are many ways to plot correlograms, the simplest one being is
    giving a data frame argument to the standard `plot` function.

    `pairs` is extremely flexible, because it allows you to specify your
    plotting function. In the above, we use the `smoothScatter` to produce
    plots in the upper triangle, and a custom text printing function in the
    lower panel to show the correlation and p-value.

 4. How many genes, in each contrast, have absolute log2 fold change
    greater than 1 and adjusted p values smaller than 0.01? Write a
    function that calculates this number from a list of contrasts given two
    paramters, lfc.thr and qval.thr.


    **Solution:**

    ```
    cid <- 1
    sel.contrast <- abs(files[[cid]]$log2FoldChange) > 1
    sel.pval     <- files[[cid]]$padj < .01
    sum(sel.contrast & sel.pval, na.rm=TRUE)

    ## this function takes a data frame
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
    ```

    **Explanations:**

    Remember that `sum(x)` returns the number of TRUE values in a logical
    vector (because TRUE == 1 and FALSE == 0). Remember that we need to
    take care of the NA's with `na.rm=TRUE`.

    We define two functions in the above: the first one does the actual
    work, and takes a data frame. The second is a wrapper around `sapply`
    and not really necessary, but I put it here because that was the task
    given. However, examine how default values are defined in the function
    definitions and how they are passed from one function to another.
