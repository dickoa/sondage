test_that("expected_hits errors on WOR design", {
  s <- equal_prob_wor(10, 3)
  expect_error(expected_hits(s), "with-replacement")
})

test_that("joint_inclusion_prob errors on WR design", {
  s <- equal_prob_wr(10, 3)
  expect_error(joint_inclusion_prob(s), "without-replacement")
})

test_that("joint_inclusion_prob.default errors on non-design object", {
  expect_error(joint_inclusion_prob(list()), "WOR design")
  expect_error(joint_inclusion_prob(1:5), "WOR design")
})

test_that("joint_expected_hits errors on WOR design", {
  s <- equal_prob_wor(10, 3)
  expect_error(joint_expected_hits(s), "with-replacement")
})

test_that("joint_expected_hits.default errors on non-design object", {
  expect_error(joint_expected_hits(list()), "WR design")
})

test_that("sampling_cov.default errors on non-design object", {
  expect_error(sampling_cov(list()), "sampling functions")
})
