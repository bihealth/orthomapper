## code to prepare `taxonomy`  table for the orthology DB goes here
## we query NCBI taxonomy database to get information for all the species
## in the orthology database


require(orthomapper)
sp <- species_all()
sptxt <- paste(sp, collapse=",")

require(xml2)
url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=%s&retmode=xml"
foo <- read_xml(sprintf(url, sptxt))

xmlsp <- xml_children(foo)

scinames <- sapply(xmlsp, function(x) as_list(x)$ScientificName[[1]])

commonnames <- sapply(xmlsp, function(x) {
  x <- as_list(x)$OtherNames

  names <- c(unlist(x$CommonName), unlist(x$GenbankCommonName))
  paste(names, collapse="; ")
})

lineages <- sapply(xmlsp, function(x) as_list(x)$Lineage[[1]])

species <- data.frame(TaxID=sp, name=scinames, common_names=commonnames, lineage=lineages, stringsAsFactors=FALSE)

output_file <- "inst/extdata/orthologDB.sqlite"
.con <- DBI::dbConnect(RSQLite::SQLite(), dbname=output_file)
DBI::dbWriteTable(.con, name="species", value=species)
DBI::dbDisconnect(.con)




