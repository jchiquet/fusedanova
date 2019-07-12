Univariate hierarchical clustering at large scale
================

## Description

Fast implementations (in O(n log(n))) of two univariate hierarchical
clustering methods (1 dimensional Ward hierarchical clustering and l1
clusterpath). Resulting hierarchies are exported as hclust object.

## Installation

``` r
devtools::install_github("jchiquet/univarclust").
```

## References

Chiquet J, Gutierrez P, Rigaill G: *Fast tree inference with weighted
fusion penalties*, **Journal of Computational and Graphical Statistics**
205–216, 2017. [PDF
version](http://www.tandfonline.com/doi/abs/10.1080/10618600.2015.1096789?journalCode=ucgs20)

Toby Hocking, Jean-Philippe Vert, Francis Bach, and Armand Joulin.
*Clusterpath: an Algorithm for Clustering using Convex Fusion
Penalties*, Proceedings of the 28th International Conference on Machine
Learning (ICML-11), 2011.
