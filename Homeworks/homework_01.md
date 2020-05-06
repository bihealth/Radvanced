# Homework 1

 1. Use `read.table()` to load the counts data from the file
    [Data/covid_data/counts.tsv](Data/covid_data/counts.tsv). How many
    genes and columns are there? How many samples? What are the sample
    names? What data type is the result? What are the gene IDs?

 2. Load the file `annotation_all.tsv`. What is it? How can you use it?

 3. Plot density of the sample 1 using functions `plot` and `density`. How
    can you add density plots for other samples with the `lines` function?
    Can you automatize it with `lapply`? How does the density plot look
    like and why? What happens when you add the parameter `log="x"` to the
    `plot()` function?

 4. Load the file `DESeq2.all.rld.model.csv`. Plot the first sample from
    this file against the first sample from the counts file. What do you
    notice?

 5. Go to the GEO and find the data set GSE147507. Read the description and
    find the expression matrix. Download the expression matrix and try to
    load it into R (this will not be straightforward! please inspect the
    file carefully).

 6. Based on the GEO information and the sample names from the matrix from
    step 1, create a data frame with sample names, cell type, treatment
    (virus name or "mock").

 7. With the help of the annotation, locate the rows corresponding to genes
    CXCL5, IL6 and IFNG in the matrices from steps 1 and 4.

 7. Using `boxplot()` and (independently) `beeswarm`, plot the expression
    of these genes in the different experimental groups. On the x axis
    there should be 8 groups (NHBE mock, NHBE Sars Cov2, A549 mock, A549
    Sars Cov2... etc), and on the y axis the log counts.
