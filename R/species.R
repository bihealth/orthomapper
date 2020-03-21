#' Return the species DBI table
#'
#' Return the species DBI table
#'
#' This file contains information about annotationDBI org.*db files of
#' which orthomapper is aware. For the taxon IDs in this database, the
#' respective org.*db will be loaded automatically if necessary.
#' @export
speciesDBITable <- function() {

  return(dbifiles)

}


#' List taxon IDs in the orthomapper DB
#' 
#' List taxon IDs in the orthomapper DB
#' @return a numeric vector of taxon IDs
species_all <- function() {

  fmt <- "SELECT DISTINCT TaxID1 FROM orthologs"
  map1 <- DBI::dbGetQuery(con, fmt)
  fmt <- "SELECT DISTINCT TaxID2 FROM orthologs"
  map2 <- DBI::dbGetQuery(con, fmt)

  return(sort(unique(unlist(c(map1, map2)))))
}
