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