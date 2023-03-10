% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orthomap.R
\name{orthologs_all}
\alias{orthologs_all}
\title{Find all orthologs to a gene}
\usage{
orthologs_all(gene)
}
\arguments{
\item{gene}{a vector of gene IDs}
}
\value{
list with each element corresponding to one gene from the `gene`
        vector. Each element is a data frame with a taxon ID and gene ID column.
}
\description{
Find all orthologs across all species to a given gene
}
\examples{
\dontrun{
orthologs_all(52024)
}
}
