Class 4: R Basics
================

In-Class exercise as a part of BGGN213 W19 at UCSD Instructor: Barry Grant

Exercise 2: Simple Calculations

Starting with attempting basic arthmetic

Can do addition:

``` r
5 + 3
```

    ## [1] 8

Subtraction:

``` r
19  -4
```

    ## [1] 15

Multiplication

``` r
4 * 7 
```

    ## [1] 28

Division

``` r
15 / 3
```

    ## [1] 5

Exercise 3: Saving your answers using object assignmnet

Assign the result of a calculation to an object. View object.

``` r
x <- 2 * 7 * 9
x
```

    ## [1] 126

``` r
# NOTE use <- and not = which will cause issues later
```

Exercise 4: Calling Functions

R has a lot of built in functions that can make calculations quick and easy. Here are a few examples.

``` r
seq(1,10)
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

If you don't know how a function works then you can use help(function) to get information about it

``` r
help(seq)
```

Using additional arguments to seq

``` r
seq(1,10, by =2)
```

    ## [1] 1 3 5 7 9

Exercise 5: How to get help in R

``` r
help(log)
?log
```

If you don't know which function you need you can do the following

``` r
help.search("cross tabulate")
??"cross tabulate"
```

Each function has an example file

``` r
example(log)
```

    ## 
    ## log> log(exp(3))
    ## [1] 3
    ## 
    ## log> log10(1e7) # = 7
    ## [1] 7
    ## 
    ## log> x <- 10^-(1+2*1:9)
    ## 
    ## log> cbind(x, log(1+x), log1p(x), exp(x)-1, expm1(x))
    ##           x                                                    
    ##  [1,] 1e-03 9.995003e-04 9.995003e-04 1.000500e-03 1.000500e-03
    ##  [2,] 1e-05 9.999950e-06 9.999950e-06 1.000005e-05 1.000005e-05
    ##  [3,] 1e-07 1.000000e-07 1.000000e-07 1.000000e-07 1.000000e-07
    ##  [4,] 1e-09 1.000000e-09 1.000000e-09 1.000000e-09 1.000000e-09
    ##  [5,] 1e-11 1.000000e-11 1.000000e-11 1.000000e-11 1.000000e-11
    ##  [6,] 1e-13 9.992007e-14 1.000000e-13 9.992007e-14 1.000000e-13
    ##  [7,] 1e-15 1.110223e-15 1.000000e-15 1.110223e-15 1.000000e-15
    ##  [8,] 1e-17 0.000000e+00 1.000000e-17 0.000000e+00 1.000000e-17
    ##  [9,] 1e-19 0.000000e+00 1.000000e-19 0.000000e+00 1.000000e-19

Excercise 6: Vectors and Indexing

Vectors are a means for storing data in R

``` r
length(3.1)
```

    ## [1] 1

To combine values into vectors use c()

``` r
x <- c(121, 4.56, 8.235)
x
```

    ## [1] 121.000   4.560   8.235

``` r
y <- c(898, 1.543, 705)
y
```

    ## [1] 898.000   1.543 705.000

Vectorization means that R loops over the values in vectors. Allows us to use arithmetic functions to vectors of the same length.

``` r
x + y
```

    ## [1] 1019.000    6.103  713.235

``` r
x - y 
```

    ## [1] -777.000    3.017 -696.765

``` r
x / y 
```

    ## [1] 0.13474388 2.95528192 0.01168085

``` r
sqrt(x)
```

    ## [1] 11.000000  2.135416  2.869669

``` r
round(sqrt(x), 3)
```

    ## [1] 11.000  2.135  2.870

``` r
log(x)/2 + 1
```

    ## [1] 3.397895 1.758661 2.054197

Vector Indexing: We can use indexing to access specific values within a vector.

``` r
x <- c(56, 95.3, 0.4)
x[2]
```

    ## [1] 95.3

Index 1 responds to the first value in a vector

``` r
x[1]
```

    ## [1] 56

Note on reproducibility

Use sessionInfo() to record the version of R and the packages used to make sure everything is reproducible

``` r
sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.4
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.2  magrittr_1.5    tools_3.5.2     htmltools_0.3.6
    ##  [5] yaml_2.2.0      Rcpp_1.0.1      stringi_1.4.3   rmarkdown_1.12 
    ##  [9] knitr_1.22      stringr_1.4.0   xfun_0.6        digest_0.6.18  
    ## [13] evaluate_0.13
