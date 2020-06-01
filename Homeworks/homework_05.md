 1. Download the data set `WIID_19Dec2018.xlsx`. This is a spreadsheet
    containing worlid income inequality data. 

 1. Load the spreadsheet in R using the `readxl` package. How can you find
    out how many sheets are in the file (without opening it in Excel, of
    course!). Hint: take a look at what functions are there in the `readxl`
    package.
 
 2. Examine the tibble produced by the `read_excel` function (why
    `read_excel` and not `read_xlsx`? Read the manual). `q1..q5` are the
    quantiles and `d1..d10` are the deciles. I.e., `q1` equal to 9 means
    that the lower income 20% of the population owns 9% of the overall
    wealth. The Gini coefficient ranges from 0 (perfect equality) to 100%
    (perfect inequality).

 3. Which countries have the highest / lowest Gini coefficient? Which
    countries had the highest / lowest Gini coefficient in 2016? Where is
    Germany on that scale?

 2. Construct a tidyverse pipe with the following steps:

     1. Select only entries from 2016
     2. Remove duplicate entries for each country (you can use the
        `duplicated` function for that)
     3. Select the columns id, country and the quantile columns
     4. Use the `gather` function to create a long version of the table 
     5. Remove rows containing NA's using the `drop_na` function

 3. Using `group_by`, `summary` and `arrange` find out for which year there
    are the most data? For which year the overall average gini coefficient
    is highest / lowest? (mind the NA's!). Plot the average Gini
    coefficient over time.
