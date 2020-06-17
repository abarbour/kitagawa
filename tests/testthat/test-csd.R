skip('no cross spectrum yet -- until psd updated etc')
context('Cross spectrum')

set.seed(1222)
x <- xna <- rnorm(100)
xna[sample(100,10)] <- NA

test_that("NA values stop execution",{
  expect_error(cross_spectrum(cbind(x,xna)))
  expect_error(cross_spectrum(cbind(xna,x)))
  expect_error(cross_spectrum(cbind(xna,xna)))
})
test_that("Types and classes: MT vs WOSA",{
  XY <- cbind(x,x^2)
  expect_type(cross_spectrum(XY), 'list')
  expect_type(cross_spectrum(XY, k=NULL), 'list')
  expect_is(cross_spectrum(XY), 'mtcsd')
  expect_is(cross_spectrum(XY, k=NULL), 'wosacsd')
})