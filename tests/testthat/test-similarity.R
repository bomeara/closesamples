test_that("get similarity works", {
  phy <- ape::rcoal(10)
  focal_species <- phy$tip.label[1:4]
  similarity_matrix <- GetSimilarity(phy, focal_species)
  expect_equal(nrow(similarity_matrix), 10)
  expect_equal(ncol(similarity_matrix), 4)
})

test_that("get data.frame of closest taxa works", {
  phy <- ape::rcoal(10)
  focal_species <- phy$tip.label[1:5]
  closest.df <- GetClosestSamples(n=3, phy, focal_species)
  expect_equal(nrow(closest.df), 3)
})

test_that("replacing works", {
  ntax=10
  phy <- ape::rcoal(ntax)
  focal_species <- phy$tip.label
  closest.df <- GetClosestSamples(n=ntax, phy, focal_species, replace_full=FALSE, replace_feasible=FALSE)
})
