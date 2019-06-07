Untitled
================

``` r
data$ali[,41]
```

    ##     P53_wt P53_mutant 
    ##        "D"        "L"

``` r
start.ind <- 41-8
end.ind <- 41 + 8

data$ali[, start.ind:end.ind]
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## P53_wt     "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "D"  "D"   "L"   "M"  
    ## P53_mutant "S"  "P"  "L"  "P"  "S"  "Q"  "A"  "M"  "L"  "D"   "L"   "M"  
    ##            [,13] [,14] [,15] [,16] [,17]
    ## P53_wt     "L"   "S"   "P"   "D"   "D"  
    ## P53_mutant "L"   "S"   "P"   "D"   "D"

``` r
missmat<- conserv(data, method = "identity")
missmat <- which(missmat < 1)
```

Find the missmatches that are not just a gap

``` r
gaps <- gap.inspect(data)
gap.inds <- gaps$t.inds


tumor_mut <- missmat[!missmat %in% gap.inds]
```

``` r
#data$ali[i, tumor_mut]
  
ids <- paste(data$ali[1, tumor_mut], tumor_mut, data$ali[2,tumor_mut], sep = "")
```

``` r
start.ind <- tumor_mut - 8
end.ind <- tumor_mut + 8

tumor <- NULL
for(i in 1:length(start.ind)){
  tumor <- seqbind(tumor, data$ali[2, start.ind[i]:end.ind[i]])
}
```
