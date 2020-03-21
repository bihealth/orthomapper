# Orthomaper

Orthomapper uses an NCBI database of orthologs to retrieve orthologs
between two taxons. 

# Installation

``` bash
git clone git@cubi-gitlab.bihealth.org:january.weiner/orthomapper.git
R CMD INSTALL orthomapper
```

Note: the package is quite large (~130 MB), as it contains the orthology DB.

## Example

Here are some examples of things you can do with orthomapper:

``` r
## get all orthologs between mouse and rat
mouserat <- orthologs(10090, 10116)

## get the SYMBOL for gene 52024, infer taxon automatically
entrez_annotate(52024)

## get more information on the gene ID 52024
entrez_annotate(52024, taxon=10090, c("SYMBOL", "GENENAME"))

## get all orthologs of the gene ID 52024
ankrd22.o <- orthologs_all(52024)
```

For more information, refer to the package vignette.
