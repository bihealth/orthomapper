#' Create orthology database from flat file
#'
#' Convert orthology database stored as a tab separated values file to SQLite
#' database
#'
#' This function was used to derive the sqlite database of orthologs
#' included in the orthomapper package. The original source for the orthologs
#' was the file gene_orthologs.gz downloaded from
#' the NCBI gene FTP site, https://ftp.ncbi.nlm.nih.gov/gene/DATA/.
#' @param input_file Path to input file (flat text or gzipped)
#' @param output_file Path to SQLite output file (must not exist!)
#' @param fields Columns in input_file corresponding to source taxon ID, source gene ID, target taxon ID and target gene ID
#' @param info named list holding additional information (comment) which should be stored in the database
#' @import DBI RSQLite
#' @importFrom utils read.table
#' @examples
#' \dontrun{
#' 
#' }
#' @export
orthologyDBFromFile <- function(input_file, output_file, fields=c(1,2,4,5), info=NULL) {

  if(file.exists(output_file))
    stop(sprintf("Cowardly refusing to overwrite existing file %s", output_file))

	orthologs <- read.table(input_file, header=T, comment.char="", sep="\t")
  orthologs <- orthologs[ , fields ]
  colnames(orthologs) <- c("TaxID1", "GeneID1", "TaxID2", "GeneID2")

  if(!is.null(info)) {
    if(is.null(names(info))) stop("`info` must be a named list")
  }

  .con <- DBI::dbConnect(RSQLite::SQLite(), dbname=output_file)
  DBI::dbWriteTable(.con, name="orthologs", value=orthologs)

  meta <- data.frame(ID=c("date", "source"), Vals=c(as.character(Sys.Date()), input_file), stringsAsFactors=FALSE)

  if(!is.null(info)) {
    .info <- simplify2array(info)
    meta <- rbind(meta, data.frame(ID=names(.info), Vals=as.character(.info)), stringsAsFactors=FALSE)
  }

  DBI::dbWriteTable(.con, name="meta", value=meta)
  DBI::dbDisconnect(.con)
}

