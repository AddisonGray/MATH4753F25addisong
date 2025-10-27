test_that("myncurve returns correct mu", {
  expect_equal(myncurve(10,5,6)$mu, 10)
})
test_that("myncurve returns correct sigma", {
  expect_equal(myncurve(10,5,6)$sigma, 5)
})
test_that("myncurve returns correct area", {
  expect_equal(round(myncurve(10,5,6)$area, 3), round(pnorm(6,10,5), 3))
})
