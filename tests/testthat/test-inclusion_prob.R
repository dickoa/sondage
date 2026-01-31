test_that("inclusion_prob returns correct length", {
  pik <- inclusion_prob(c(10, 20, 30, 40), n = 2)
  expect_length(pik, 4)
})

test_that("inclusion_prob sums to n", {
  pik <- inclusion_prob(c(10, 20, 30, 40), n = 2)
  expect_equal(sum(pik), 2)

  pik <- inclusion_prob(1:100, n = 25)
  expect_equal(sum(pik), 25)
})

test_that("inclusion_prob values are proportional to size", {
  size <- c(10, 20, 30, 40)
  pik <- inclusion_prob(size, n = 2)

  expect_equal(pik[2] / pik[1], 2, tolerance = 1e-6)
  expect_equal(pik[3] / pik[1], 3, tolerance = 1e-6)
  expect_equal(pik[4] / pik[1], 4, tolerance = 1e-6)
})

test_that("inclusion_prob handles certainty selections", {
  # One very large unit
  size <- c(1, 1, 1, 100)
  pik <- inclusion_prob(size, n = 2)

  expect_equal(pik[4], 1) # Large unit gets probability 1
  expect_equal(sum(pik[1:3]), 1) # Remaining units share 1
})

test_that("inclusion_prob handles multiple certainty selections", {
  size <- c(1, 100, 200, 300) # Three large units
  pik <- inclusion_prob(size, n = 3)

  expect_true(pik[2] == 1 || pik[3] == 1 || pik[4] == 1)
})

test_that("inclusion_prob returns values in [0, 1]", {
  set.seed(42)
  size <- runif(100, 1, 1000)
  pik <- inclusion_prob(size, n = 30)

  expect_true(all(pik >= 0))
  expect_true(all(pik <= 1))
})

test_that("inclusion_prob handles equal sizes", {
  size <- rep(1, 10)
  pik <- inclusion_prob(size, n = 5)

  expect_equal(pik, rep(0.5, 10))
})

test_that("inclusion_prob warns about negative values", {
  expect_warning(
    pik <- inclusion_prob(c(-1, 2, 3, 4), n = 2),
    "negative"
  )
})

test_that("inclusion_prob rejects invalid input", {
  expect_error(inclusion_prob(1:10, n = 15), "cannot exceed")
  expect_error(inclusion_prob(1:10, n = -1), "non-negative")
  expect_error(inclusion_prob(1:10, n = NA), "not NA")
})

test_that("inclusion_prob works with up_maxent", {
  size <- c(500, 1200, 800, 3000, 600)
  pik <- inclusion_prob(size, n = 3)

  set.seed(42)
  idx <- up_maxent(pik)
  expect_length(idx, 3)
})
