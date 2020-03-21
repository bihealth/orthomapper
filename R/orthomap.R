con <- NULL

.onLoad <- function(libname, pkgname) {
  dbfile <- system.file("extdata", "orthologDB.sqlite", package="orthomapper", mustWork=TRUE)
  con <<- DBI::dbConnect(RSQLite::SQLite(), dbname=dbfile)
}

.onUnload <- function(libpath) {
  DBI::dbDisconnect(con)
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
#' @param ids entrez gene IDs from taxon1
#' @param taxon1 source taxon ID
#' @param taxon2 target taxon ID
#' @return A data frame with two columns
#' @import DBI
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


## below is the old implementation not using self join. Slightly faster,
## but less elegant.

#' All ortholog pairs between two taxons
#'
#' All ortholog pairs between two taxons
#' @param ids entrez gene IDs from taxon1
#' @param taxon1 source taxon ID
#' @param taxon2 target taxon ID
#' @return A data frame with two columns
#' @import DBI
.orthologs <- function(taxon1, taxon2) {
  
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


.getdbi <- function(taxon) {
  m <- match(taxon, dbifiles$taxonID)
  if(is.na(m)) {
    stop(sprintf("Can't find a matching org.db for taxon %d automatically", taxon))
  }

  orgdb <- dbifiles$DB[m]
  return(orgdb)
}

.loaddbi <- function(orgdb) {
  if(!require(orgdb, character.only=TRUE, quietly=TRUE)) {
    stop(sprintf("Cannot load package %s, is it installed?", orgdb))
  }

  dbi <- eval(parse(text=orgdb))
  return(dbi)
}

#' Return annotation information for given Entrez IDs
#'
#' Return annotation information for given Entrez IDs
#' @param genes a numeric vector of Entrez gene IDs
#' @param taxon a single value for a taxon. If NULL, orthomapper will
#'              attempt to identify the taxon (it is much faster if you specify it)
#' @param column a character vector with the columns requested from the annotation database
#' @param orgdb name of the annotation package; inferred automatically for
#'              species present in the output of the `speciesDBITable()` function.
#' @return a data frame with `length(column)+1` columns; first column is the entrez ID
#' @examples
#' ## get the SYMBOL for gene 52024, infer taxon automatically
#' entrez_annotate(52024)
#' 
#' ## get more information on the gene ID 52024
#' entrez_annotate(52024, taxon=10090, c("SYMBOL", "GENENAME"))
#' @export
entrez_annotate <- function(genes, taxon=NULL, column="SYMBOL", orgdb=NULL) {

  if(is.null(taxon)) 
    taxon <- find_taxon(genes[1])[1,2]

  if(is.null(taxon) || is.na(taxon))
    stop("Cannot identify taxon!")

  if(is.null(orgdb)) orgdb <- .getdbi(taxon)
  dbi <- .loaddbi(orgdb)

  cols <- columns(dbi)

  if(any(missing <- ! column %in% cols)) {
    stop(sprintf("Following columns not defined in %s: %s",
      orgdb, paste(column[missing], collapse=", ")))
  }

  genes2 <- as.character(genes)
  ret <- lapply(column, function(cc) {
    mapIds(dbi, genes2, column=cc, keytype="ENTREZID", multiVals="first")
  })
  ret <- Reduce(cbind, ret)

  ret <- data.frame(genes, ret, stringsAsFactors=FALSE)
  colnames(ret) <- c("Entrez", column)

  return(ret)
}

#' Search for symbols similar to a given string
#'
#' Search for symbols similar to a given string
#'
#' Basically, uppercase the requested symbol, remove all digits at the end
#' and see whether anything matches.
#' @param orgdb name of the annotation package; inferred automatically for
#'              species present in the output of the `speciesDBITable()` function.
#' @param taxon a single value for a taxon
#' @return a data frame with 2 columns; first column is the entrez ID, the second the symbol, third the gene name.
#' @examples
#' ## find a gene symbol by similarity search
#' similar_symbol("GBP5", 10090)
#' @export
similar_symbol <- function(x, taxon, orgdb=NULL) {

  if(is.null(orgdb)) orgdb <- .getdbi(taxon)
  dbi <- .loaddbi(orgdb)

  keys <- keys(dbi, 'ENTREZID')
  res <- select(dbi, keys=keys, columns=c("ENTREZID", "SYMBOL", "GENENAME"))

  x <- gsub("[0-9]*$", "", x)
  sel <- grep(toupper(x), toupper(res$SYMBOL))
  return(res[sel, ])
}
