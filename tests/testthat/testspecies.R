context("Species DB")


test_that("Orthomapper has species DB", {
  a <- species_info()
  expect_true(nrow(a) > 300)
  a <- species_info(9606)
  expect_equal(nrow(a), 1)
  expect_equal(a$TaxID, 9606)
  expect_match(a$name, "Homo sapiens")
})

test_that("Species DB matches the orthologsDB", {
  n1 <- length(species_all())
  n2 <- nrow(species_info())
  expect_equal(n1, n2)
})

test_that("Orthomapper can find species", {
  s <- species_search("human")
  expect_equal(nrow(s), 1)
})


