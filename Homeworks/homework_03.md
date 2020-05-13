
 1. The standard deviation in a population is defined as 

    <img
    src="https://render.githubusercontent.com/render/math?math=\sqrt{\Sum_i^n\frac{(x_i-\overline{x})^2}{n}">

    
    `sqrt(sum((x - mean(x))^2)/length(x))`. However, in statistics we use 
    `sqrt(sum((x - mean(x))^2)/(length(x)-1))` when calculating the
    *estimator* of standard deviation based on a sample. Using a simulation
    approach, show that indeed the latter is more suitable and that as
    `length(x)` grows, the two values become closer. Here is one possible algorithm:
    
     1. For a given N, sample the normal distribution a thousand times (with known SD as a
        parameter; use the `rnorm` function), each time drawing N numbers.
        
     2. Calculate the SD with both methods. 

     2. Plot 
    
     Which one is closer to the real SD of the population
    from which the sample was drawn?
