
 1. Load all the contrast files. Can you do it automatically? Look up the
    function `list.files()` and its parameter `pattern`. Use `lapply` to
    load all the files.

    **Solution:**

    ```
    files.list <- list.files("Data/covid_data/", pattern="contrast*", full.names=TRUE)
    files <- lapply(files.list, read.table, header=TRUE)
    ```

    **Explanations:**



 2. Create a matrix containing only the log2 fold changes from all
    contrasts.  Either use `lapply` and `Reduce`, or try `sapply` which can
    do some automatic transformations of the output.

 3. Plot the correlation between two contrasts. Can you make a plot that
    shows all the correlations in one go? Either take a look at the web
    site
    [https://www.r-graph-gallery.com/correlogram.html](https://www.r-graph-gallery.com/correlogram.html)
    or (which is harder, but more flexible) use the `pairs` or `ggpairs`
    functions (`ggpairs` is in the GGally package).

 4. How many genes, in each contrast, have absolute log2 fold change
    greater than 1 and adjusted p values smaller than 0.01? Write a
    function that calculates this number from a list of contrasts given two
    paramters, lfc.thr and qval.thr.
    
