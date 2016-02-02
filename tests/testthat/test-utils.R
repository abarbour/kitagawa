context('Utility functions')

Transmiss <- 1e-6
Storativ <- 1e-6
omeg <- 10**seq(-3,3, by=0.1) # Hz
dep <- 30 # meters

Vw <- 1
Rs <- 0.1
Ku <- 1e9
B <- 0.9

frq <- omega_norm(omeg, Storativ/Transmiss, z = dep, invert = TRUE)

sresp <- well_response(frq, Transmiss, Storativ, Vw, Rs, Ku, B)
oresp <- open_well_response(frq, Transmiss, Storativ, z=dep)

test_that("classes are correct",{
  expect_true(is.wrsp(sresp))
  expect_true(is.owrsp(oresp))
})

test_that("class checking functional",{
  expect_false(is.wrsp(data.frame()))
  expect_false(is.owrsp(data.frame()))
  
  expect_false(is.wrsp(matrix()))
  expect_false(is.owrsp(matrix()))

  expect_false(is.wrsp(list()))
  expect_false(is.owrsp(list()))
})