Class 12: Structural Bioinformatics
================

Section 1: In silico docking of drugs to HIV-1 protease
=======================================================

This exercise is a continuation of Class 11. We will be looking at the structure of HIV-1 protease.

``` r
library(bio3d)

file.name <- get.pdb("1hsg")
```

    ## Warning in get.pdb("1hsg"): ./1hsg.pdb exists. Skipping download

Next use the read.pdb()function to read this PDB file into R so we can prepare it for further analysis

``` r
hiv <-read.pdb(file.name)
```

What does this pdb file look like?

``` r
hiv
```

    ## 
    ##  Call:  read.pdb(file = file.name)
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

What are the two non-protein values in this structure? HOH (127), MK1 (1) water! and probably a ligand

We can separate all of these components by using trim.pdb()

``` r
# separate the atoms associated with protein and ligand
# save them separately 

prot <- trim.pdb(hiv, "protein")
lig <- trim.pdb(hiv, "ligand")

# write them as pdb files to my local directory

write.pdb(prot, file="1hsg_protein.pdb")
write.pdb(lig, file="1hsg_ligand.pdb")
```

Section 2: Docking ligands into HIV-1 protease
==============================================

Load in our docking results and look at them

``` r
#vina is having issues 

res <- read.pdb("all.pdbqt", multi=TRUE)

write.pdb(res, "results.pdb")
```

Calculate the root mean square distance between the docking results and the pdb file

``` r
ori <- read.pdb("ligand.pdbqt")

rmsd(ori, res)
```

    ##  [1]  0.697  4.195 11.146 10.606 10.852 10.945 10.945  3.844  5.473  4.092
    ## [11] 10.404  5.574  3.448 11.396  6.126  3.848  8.237 11.196 10.981 11.950

Section 3: Exploring confirmational dynamics
============================================

Normal mode analysis is a common simulation technique to look at biomolecules

First we will download a pdb file of hen egg white lysozyme (PDB id 1hel)

``` r
library(bio3d)

pdb <- read.pdb("1HEL")
```

    ##   Note: Accessing on-line PDB file

``` r
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.016 seconds.
    ##  Diagonalizing Hessian...    Done in 0.089 seconds.

``` r
plot(modes, sse=pdb)
```

![](Class_12-_Structural_Bioinformatics_files/figure-markdown_github/unnamed-chunk-7-1.png)

Use the function mktrj() to generate a trajectory PDB file by interpolating along a given normal mode:

``` r
mktrj(modes, mode=7, file="nma_7.pdb")
```

This can be openend in VMD to look at it.

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
    ## other attached packages:
    ## [1] bio3d_2.3-4
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] compiler_3.5.2  magrittr_1.5    parallel_3.5.2  tools_3.5.2    
    ##  [5] htmltools_0.3.6 yaml_2.2.0      Rcpp_1.0.1      stringi_1.4.3  
    ##  [9] rmarkdown_1.12  grid_3.5.2      knitr_1.22      stringr_1.4.0  
    ## [13] xfun_0.6        digest_0.6.18   evaluate_0.13
