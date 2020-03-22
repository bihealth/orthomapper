context("Mapping orthologs")

test_that("orthomap returns orthologs", {
  expect_equal(nrow(orthologs(10090, 10116)), 15806)
  expect_equal(nrow(orthomap(c(229898, 52024), 10090, 9606)), 2)
})
  

test_that("Orthologs independent of order of taxa", {

  a1 <- orthologs(10090, 10116)
  a2 <- orthologs(10116, 10090)
  expect_true(setequal(paste(a1[,1], a1[,2]), paste(a2[,2], a2[,1])))

})

test_that("orthomap works", {
  a1 <- orthomap(52024, 10090, 9606)
  expect_equal(nrow(a1), 1)
  expect_equal(a1[,2], 118932)
})

test_that("orthologs_all works", {
  a1 <- orthologs_all(52024)

  expect_true(is.list(a1))
  expect_equal(length(a1), 1)

  expect_equal(nrow(a1[[1]]), 290)
  expect_identical(colnames(a1[[1]]), c("Taxon", "Gene"))

})
