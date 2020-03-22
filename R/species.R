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


## search for a single taxon name
.get_taxon <- function(name, fuzzy=TRUE) {

  if(toupper(name) %in% toupper(dbifiles$moniker)) 
    return(dbifiles$taxonID[ match(toupper(name), toupper(dbifiles$moniker)) ])

  if(any(sel <- grepl(paste0("^", name, "$"), dbifiles$species, ignore.case=TRUE))) {
    return(dbifiles$taxonID[ which(sel)[1] ])
  }

  if(!fuzzy) return(NA)

  search <- species_search(name, columns=c("name", "common_names"))

  if(nrow(search) == 0) return(NA)
  if(nrow(search) == 1) return(search$TaxID)

  warning(sprintf("Multiple species match %s, returning the first one (%s[%s])",
    name, search$name[1], search$TaxID[1]))

  return(search$TaxID[1])
}


#' Match taxon based on a string
#'
#' Return an NCBI taxon ID based on species name 
#' or moniker (such as "Hs" for Homo sapiens or "Mm" for mouse).
#' Note that if multiple species match the strings, only the first one will
#' be returned.
#' @return named numeric ID vector of the taxons
#' @seealso species_search
#' @param name character vector with species names or monikers
#' @export
get_taxon <- function(name) {

  ret <- sapply(name, .get_taxon)
  names(ret) <- name
  return(ret)
}


#' List taxon IDs in the orthomapper DB
#' 
#' List taxon IDs in the orthomapper DB
#' @return a numeric vector of taxon IDs
#' @export
species_all <- function() {

  fmt <- "SELECT DISTINCT TaxID1 FROM orthologs"
  map1 <- DBI::dbGetQuery(con, fmt)
  fmt <- "SELECT DISTINCT TaxID2 FROM orthologs"
  map2 <- DBI::dbGetQuery(con, fmt)

  return(sort(unique(unlist(c(map1, map2)))))
}



#' Retrieve species info 
#'
#' Retrieve species info either for all species in the orthology database
#' or for the species requested throught the taxonID parameter.
#' @param taxonID numeric vector of taxonomy IDs. If NULL, retrieve all species. Otherwise, find the
#'                species selected
#' @param show.lineage if TRUE, show also the lineage column (not often needed)
#' @examples
#' species_info(c(9606, 10090))
#' @export
species_info <- function(taxonID=NULL, show.lineage=FALSE) {

  if(is.null(taxonID)) {
    query <- "SELECT * FROM species"
  } else {
    query <- "SELECT * FROM species WHERE TaxID IN (%s)"
    query <- sprintf(query, paste(taxonID, collapse=", "))
  }

  ret <- DBI::dbGetQuery(con, query)

  if(!is.null(taxonID)) {
    ret <- ret[ match(taxonID, ret$TaxID), ]
    ret$TaxID <- taxonID
  }
  rownames(ret) <- NULL

  if(!show.lineage) ret$lineage <- NULL

  return(ret)
}


#' Find all species matching a pattern
#'
#' Find all species matching a pattern
#'
#' @param pattern regular expression to serch for
#' @param columns columns from the species table to search through
#' @param ignore.case if TRUE (default) do a case-insensitive search
#' @param show.lineage if TRUE, show also the lineage column (not often needed)
#' @examples
#' ## find Danio rerio
#' species_search("danio")
#' ## find all birds
#' species_search("aves")
#' @export
species_search <- function(pattern, columns=c("name", "common_names", "lineage"), ignore.case=TRUE,
  show.lineage=FALSE) {

  sp <- species_info(show.lineage=TRUE)

  sel <- lapply(columns, function(x) grepl(pattern, sp[,x], ignore.case=ignore.case))
  sel <- Reduce(`|`, sel)

  if(!show.lineage) sp$lineage <- NULL

  sp[sel, ]

}
