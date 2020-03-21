context("General")

library(orthomapper)

test_that("Orthomapper returns orthologs", {
  expect_equal(nrow(orthologs(10090, 10116)), 15806)
  expect_equal(nrow(orthomap(c(229898, 52024), 10090, 9606)), 2)
})
  

test_that("Orthomapper returns annotations", {
    expect_match(entrez_annotate(115362)[,2], "GBP5")
})
