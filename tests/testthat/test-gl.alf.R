test_that("number of rows",{
  expect_equal(nrow(gl.alf(platypus.gl)),1000)
})

test_that("class of return",{
expect_equal(class(gl.alf(platypus.gl)),"data.frame")
})
