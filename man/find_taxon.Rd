% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orthomap.R
\name{find_taxon}
\alias{find_taxon}
\title{Return a taxon for each of the genes}
\usage{
find_taxon(genes)
}
\arguments{
\item{genes}{numeric vector of Entrez gene IDs}
}
\value{
a data frame with two columns
}
\description{
Return a taxon for each of the genes
}
\examples{
find_taxon(c(115362, 229898))
}
