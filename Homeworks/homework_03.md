
 1. The standard deviation in a population is defined as 

    <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\frac{\sum_i=1^n(x_i-\overline{x})^2}{n}}">

    Or, in R, 
    `sqrt(sum(((x - mean(x))^2))/length(x))`. However, in statistics we use 

    <img src="https://render.githubusercontent.com/render/math?math=\sqrt{\frac{\sum_i=1^n(x_i-\overline{x})^2}{n - 1}}">

    in R: `sqrt(sum(((x - mean(x))^2))/(length(x)-1))` when calculating the
    *estimator* of standard deviation based on a sample. Using a simulation
    approach, show that indeed the latter is more suitable and that as
    `length(x)` grows, the two values become closer. Here is one possible algorithm:
    
     1. For a given N, sample the normal distribution a hundred times (with known SD as a
        parameter; use the `rnorm` function), each time drawing N numbers.
        
     2. Calculate the SD with both methods for each of the hundred samples
        and calculate the average for both estimators.

     2. Plot the averages of both estimators against N. Draw a horizontal
        line (e.g. with `abline`) at the SD that was used to generate the
        samples.
    
 2. Using the approach presented in the lecture and the *full* data set,
    create the object `ds` (DESeq object) with design `~ 0 + group`.
    Create results for the following contrasts and store the results in a
    list:

     * NHBE.SC2V vs NHBE.mock
     * A549.SC2V vs A549.mock
     * A549.RSV vs A549.RSV.mock
     * A549.IAV vs A549.IAV.mock


 2. Interaction is when you have two factors which are not independent.
    For example, consider two types of cells and a treatment vs control.
    If the difference between treatment and control is *different* between
    the cell types, then there is an interaction.
    
    The table below shows example average expression of six genes in
    treatment of NHBE and A549 cells with either a mock treatment or the
    Sars-Cov-2:


    |Gene ID   |NHBE, mock|NHBE, SCov2|A549, mock|A549, SCov2|
    |----------|----------|-----------|----------|-----------|
    |AAA       |10        |1000       |10        |1000       |
    |BBB       |10        |100        |100       |1000       |
    |CCC       |0         |0          |1000      |1000       |
    |DDD       |10        |1000       |10        |10         |
    |EEE       |10        |1000       |10        |100        |
    |FFF       |1000      |500        |1000      |2000       |
    
    
     * The gene AAA changes its expression upon treatment in both types of
       cells. 
     * The gene BBB has both an effect of cell type (it is more highly
    expressed in A549 cells) and treatment (it's expression increases
    10-fold upon treatment), but the change caused by the treatment is
    identical in both cell types.  
     * The gene CCC does not react to
    treatment and is expressed only in the A549 cells. 
     * Gene DDD reacts to
    treatment only in NHBE cells (so its change *depends* on the cell
    type). 
     * Gene EEE changes expression in both cell types, albeit to
    different extent (there is a difference in differences between mock
    and treatment). 
     * Gene FFF changes in both cell types, but in different
    directions. 
    
    **Question:** which of the genes show an interaction?

 3. Simply by filtering the results of direct comparisons between SC2V and
    mock for A549 and NHBE cells, try to find genes that roughly follow
    the six patterns in the table above. It does not need to be
    significant, and it does not have to be very precise. Create dot plots
    (like with `beeswarm`) for the selected genes.


 3. To test for interaction one has to define a specific contrast. In the
    following, we would like to test the interaction between cell type and
    the treatment with Sars-Cov-2. Use the approach with numeric contrast
    vector. The interaction is given by the equation
    
    <img src="https://render.githubusercontent.com/render/math?math=(\text{NHBE}_\text{SC2V}-\text{NHBE}_\text{mock})-(\text{A549}_\text{SC2V}-\text{A549}_\text{mock})">

    Using primary school algebra, what value (-1, 0, 1) in the contrast
    vector should be associated with the coefficients given by
    `resultsNames(ds)`? 

 4. Inspect the obtained results for an interaction test. Plot some of the
    top genes. Do you see the interaction? Look up the gene examples you
    have selected before in the results tables. What do the associated
    log2FoldChanges mean?
