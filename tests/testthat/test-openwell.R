context("Open well response")

Transmiss <- 1e-6
Storativ <- 1e-6
omeg <- 10**seq(-3,3, by=0.1) # Hz
dep <- 30 # meters

frq <- omega_norm(omeg, Storativ/Transmiss, z = dep, invert = TRUE)
resp <- open_well_response(frq, Transmiss, Storativ, z=dep)

test_that("class is correct",{
  expect_is(resp,'owrsp')
})

test_that("method dispatch appropriate",{
  
  expect_is(as.data.frame(resp),'data.frame')
  
  expect_is(print(resp),'owrsp')
  expect_message(print(resp))
  
  expect_is(summary(resp),'summary.owrsp')
  expect_is(print(summary(resp)),'summary.owrsp')
  expect_message(print(summary(resp)))
  
})