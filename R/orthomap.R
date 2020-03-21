con <- NULL

.onLoad <- function(libname, pkgname) {
  dbfile <- system.file("extdata", "orthologDB.sqlite", package="orthomapper", mustWork=TRUE)
  con <<- DBI::dbConnect(RSQLite::SQLite(), dbname=dbfile)
}

.onUnload <- function(libpath) {
  DBI::dbDisconnect(con)
}

#' Connect to an alternative orthology DB
#'
#' Use an alternative orthology DB SQLite file, for example if you want a newer version.
#' The database must have a table "orthologs" with the defined columns
#' "TaxID1", "GeneID1", "TaxID2", "GeneID2".
#' @param file path to the SQLite file with the DB
#' @return invisibly returns the DBI connection to the database
#' @seealso [orthologyDBFromFile()] on how to generate orthology DB files.
#' @export
use_orthoDB <- function(file) {

	DBI::dbDisconnect(con)
	con <<- DBI::dbConnect(RSQLite::SQLite(), dbname=file)
	invisible(con)
}

.gene_taxon <- function(g) {

  fmt <- "SELECT DISTINCT TaxID1 AS TaxID FROM orthologs WHERE GeneID1 == %s"
  map <- DBI::dbGetQuery(con, sprintf(fmt, g))

  fmt <- "SELECT DISTINCT TaxID2 AS TaxID FROM orthologs WHERE GeneID2 == %s"
  map2 <- DBI::dbGetQuery(con, sprintf(fmt, g))

  return(rbind(map, map2)[1,])
}




#' Return a taxon for each of the genes
#'
#' Return a taxon for each of the genes
#' @param genes numeric vector of Entrez gene IDs
#' @return a data frame with two columns
#' @examples
#' find_taxon(c(115362, 229898))
#' @export
find_taxon <- function(genes) {

  taxons <- lapply(genes, .gene_taxon)  
  taxons <- Reduce(rbind, taxons)
  return(cbind(data.frame(GeneID=genes, stringsAsFactors=FALSE), taxons))
}

#' All ortholog pairs between two taxons
#'
#' All ortholog pairs between two taxons
#' @param taxon1 source taxon ID
#' @param taxon2 target taxon ID
#' @return A data frame with two columns
#' @import DBI RSQLite
#' @examples
#' \dontrun{
#' mouserat <- orthologs(10090, 10116)
#'}
#' @export
orthologs <- function(taxon1, taxon2) {

  cnames <- c(sprintf("T.%d", taxon1), sprintf("T.%d", taxon2))
  
  fmt <- "SELECT GeneID1, GeneID2 FROM orthologs WHERE (TaxID1 == %s AND TaxID2 == %s)"
  query <- sprintf(fmt, taxon1, taxon2)
  map <- DBI::dbGetQuery(con, query)
  colnames(map) <- cnames

  fmt <- "SELECT GeneID2, GeneID1 FROM orthologs WHERE (TaxID2 == %s AND TaxID1 == %s)"
  query <- sprintf(fmt, taxon1, taxon2)
  map2 <- DBI::dbGetQuery(con, query)
  colnames(map2) <- cnames
  map <- rbind(map, map2)

  ## now we need to find links which go over species in the left column
  ## we use internal join for that

  fmt <- "SELECT DISTINCT o1.GeneID2, o2.GeneID2 FROM orthologs o1 JOIN
          orthologs o2 ON o1.TaxID2 == %s AND o2.TaxID2 == %s AND o1.GeneID1 == o2.GeneID1"
  query <- sprintf(fmt, taxon1, taxon2)
  map3 <- DBI::dbGetQuery(con, query)
  colnames(map3) <- cnames
       
  map <- rbind(map, map3)

  return(map)
}



.find_one_ortholog <- function(gene) {

  fmt <- "SELECT TaxID1, GeneID1 FROM orthologs WHERE GeneID2 == %s"
  query <- sprintf(fmt, gene)
  map <- DBI::dbGetQuery(con, query)
  colnames(map) <- c("Taxon", "Gene")

  fmt <- "SELECT TaxID2, GeneID2 FROM orthologs WHERE GeneID1 == %s"
  query <- sprintf(fmt, gene)
  map2 <- DBI::dbGetQuery(con, query)
  colnames(map2) <- c("Taxon", "Gene")

  ## we also need to find intermediate pairs
  fmt <- "SELECT DISTINCT o2.TaxID2, o2.GeneID2 FROM orthologs o1 JOIN orthologs o2
          ON o1.GeneID2 == %s AND o1.GeneID1 == o2.GeneID1"
  query <- sprintf(fmt, gene)
  map3 <- DBI::dbGetQuery(con, query)
  colnames(map3) <- c("Taxon", "Gene")
  return(map3)

  return(rbind(map, map2))
}



#' Find all orthologs to a gene 
#'
#' Find all orthologs across all species to a given gene
#' @param gene a vector of gene IDs
#' @return list with each element corresponding to one gene from the `gene`
#'         vector. Each element is a data frame with a taxon ID and gene ID column.
#' @examples
#' \dontrun{
#' orthologs_all(52024)
#'}
#' @export
orthologs_all <- function(gene) {
  
  ret <- lapply(gene, .find_one_ortholog)
  names(ret) <- gene
  return(ret)
}


#' Map a set of gene IDs to orthologs from another taxon
#'
#' Map a set of gene IDs to orthologs from another taxon
#' @param genes a vector with entrez IDs of genes
#' @param taxon1 source taxon ID
#' @param taxon2 target taxon ID
#' @return A data frame with two columns
#' @examples
#' orthomap(c(229898, 52024), 10090, 9606)
#' @export
orthomap <- function(genes, taxon1, taxon2) {

	map <- orthologs(taxon1, taxon2)
	match <- match(genes, map[,1])
  ret <- data.frame(genes, map[match,2])
  colnames(ret) <- c(sprintf("T.%d", taxon1), sprintf("T.%d", taxon2))
  return(ret)
}



