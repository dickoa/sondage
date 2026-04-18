test_that("print handles long single-sample vectors", {
  set.seed(1)
  s <- equal_prob_wor(100, 25)
  expect_output(print(s), "Equal prob WOR \\[srs\\]")
  expect_output(print(s), "\\.\\.\\.")
})

test_that("print handles matrix replicates with long n", {
  set.seed(1)
  s <- equal_prob_wor(100, 25, nrep = 3)
  out <- capture.output(print(s))
  expect_true(any(grepl("3 replicates", out)))
  expect_true(any(grepl("rep 1", out)))
  expect_true(any(grepl("\\.\\.\\.", out)))
})

test_that("print handles list-shaped samples (random-size batch)", {
  set.seed(1)
  s <- equal_prob_wor(20, 5, method = "bernoulli", nrep = 2)
  expect_true(is.list(s$sample))
  out <- capture.output(print(s))
  expect_true(any(grepl("2 replicates", out)))
  expect_true(any(grepl("rep 1 \\(n=", out)))
})

test_that("print list branch truncates long replicate", {
  pik <- rep(0.9, 30)
  set.seed(1)
  s <- unequal_prob_wor(pik, method = "poisson", nrep = 2)
  out <- capture.output(print(s))
  expect_true(any(grepl("\\.\\.\\.", out)))
})

test_that("print equal_prob WR matrix branch", {
  set.seed(1)
  s <- equal_prob_wr(10, 3, nrep = 2)
  out <- capture.output(print(s))
  expect_true(any(grepl("Equal prob WR", out)))
  expect_true(any(grepl("2 replicates", out)))
})
