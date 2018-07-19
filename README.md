Fused-ANOVA: a package to fit weighted fusion penalties at large scale
================

Description
-----------

Fused-ANOVA is a penalized method that solves the one-way ANOVA problem by collapsing the coefficients of K conditions. It reconstructs a balanced tree structure between the condition with a homotopy algorithm.

For a class of weights implemented here, our homotopy algorithm is in K log(K). These weights induce a balanced tree structure and simplify the interpretation of the results. The package contains an illustrating phenotypic data set: given a trait, we reconstruct a balanced tree structure and assess its agreement with the known phylogeny. More in the vignette.

Installation
------------

``` r
devtools::install_github("jchiquet/fusedanova").
```

Reference
---------

Chiquet J, Gutierrez P, Rigaill G: *Fast tree inference with weighted fusion penalties*, **Journal of Computational and Graphical Statistics** 205–216, 2017. [PDF version](http://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1096789?journalCode=ucgs20)
