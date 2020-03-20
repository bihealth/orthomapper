con <- NULL

.onLoad <- function(libname, pkgname) {
  dbfile <- system.file("extdata", "orthologDB.sqlite", package="orthomapper", mustWork=TRUE)
  con <<- DBI::dbConnect(RSQLite::SQLite(), dbname=dbfile)
}

.onUnload <- function(libpath) {
  DBI::dbDisconnect(con)
}

#' All ortholog pairs between two taxons
#'
#' All ortholog pairs between two taxons
#' @param ids entrez gene IDs from taxon1
#' @param taxon1 source taxon ID
#' @param taxon2 target taxon ID
#' @return A data frame with two columns
#' @import DBI
#' @export
orthologs <- function(taxon1, taxon2) {
  
  fmt <- "SELECT GeneID1, GeneID2 FROM orthologs WHERE (TaxID1 == %s AND TaxID2 == %s)"
  query <- sprintf(fmt, taxon1, taxon2)
  map <- DBI::dbGetQuery(con, query)
  colnames(map) <- c("Source", "Target")

  fmt <- "SELECT GeneID2, GeneID1 FROM orthologs WHERE (TaxID2 == %s AND TaxID1 == %s)"
  query <- sprintf(fmt, taxon1, taxon2)
  map2 <- DBI::dbGetQuery(con, query)
  colnames(map2) <- c("Source", "Target")
  map <- rbind(map, map2)

  ## now we need to find links which go over species in the left column

  fmt <- "SELECT GeneID1, GeneID2 FROM orthologs WHERE TaxID2 == %s"
  query <- sprintf(fmt, taxon1)
  tax1 <- DBI::dbGetQuery(con, query)
  query <- sprintf(fmt, taxon2)
  tax2 <- DBI::dbGetQuery(con, query)
  
  common <- intersect(tax1[,1], tax2[,1])

  map3 <- data.frame(Source=tax1[ match(common, tax1[,1]), 2], 
                     Target=tax2[ match(common, tax2[,1]), 2])

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
  return(rbind(map, map2))
}

#' Find all orthologs to a gene 
#'
#' Find all orthologs across all species to a given gene
#' @param gene a vector of gene IDs
#' @return list with each element corresponding to one gene from the `gene`
#'         vector. Each element is a data frame with a taxon ID and gene ID column.
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

orthomap <- function(genes, taxon1, taxon2) {

	map <- orthologs(taxon1, taxon2)
	match <- match(genes, map[,1])

  return(data.frame(Source=genes, Target=map[match,2]))
}


#' @export
entrez2symbol <- function(genes, taxon, column="SYMBOL", orgdb=NULL) {

  if(is.null(orgdb)) {
    m <- match(taxon, dbifiles$taxonID)
    if(is.na(m)) {
      stop(sprintf("Can't find a matching org.db for taxon %d automatically", taxon))
    }

    orgdb <- dbifiles$DB[m]
  }

  if(!require(orgdb, character.only=TRUE)) {
    stop(sprintf("Cannot load package %s, is it installed?", orgdb))
  }

  dbi <- eval(parse(text=orgdb))
  symbols <- mapIds(dbi, genes, column="SYMBOL", keytype="ENTREZID", multiVals="first")

  return(symbols)
}
