context("Annotations")

test_that("Orthomapper returns annotations", {
    expect_match(entrez_annotate(115362)[,2], "GBP5")

    a <- entrez_annotate(52024, taxon=10090, c("SYMBOL", "GENENAME"))
    expect_match(a$SYMBOL, "Ankrd22")
    expect_match(a$GENENAME, "ankyrin repeat domain 22")

    a <- similar_symbol("GBP5", 9606)
    expect_true(all(c("GBP1", "GBP2", "FCGBP") %in% a$SYMBOL))
    expect_true(all(c("ENTREZID", "SYMBOL", "GENENAME") %in% colnames(a)))

})

test_that("entrez_annotate handles NAs in input", {

    a <- entrez_annotate(c(52024, NA), taxon=10090, "SYMBOL")
    expect_true(is.data.frame(a))
    expect_identical(nrow(a), 2)

})
