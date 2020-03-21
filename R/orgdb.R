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
#' @importFrom AnnotationDbi columns mapIds
#' @examples
#' \dontrun{
#' ## get the SYMBOL for gene 52024, infer taxon automatically
#' entrez_annotate(52024)
#' 
#' ## get more information on the gene ID 52024
#' entrez_annotate(52024, taxon=10090, c("SYMBOL", "GENENAME"))
#' }
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
#' @param x a single string (symbol to match)
#' @param orgdb name of the annotation package; inferred automatically for
#'              species present in the output of the `speciesDBITable()` function.
#' @param taxon a single value for a taxon
#' @return a data frame with 2 columns; first column is the entrez ID, the second the symbol, third the gene name.
#' @importFrom AnnotationDbi select keys
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