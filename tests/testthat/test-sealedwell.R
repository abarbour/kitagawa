context("Sealed well response")

Transmiss <- 1e-6
Storativ <- 1e-6
omeg <- 10**seq(-3,3, by=0.1) # Hz
dep <- 30 # meters

Vw <- 1
Rs <- 0.1
Ku <- 1e9
B <- 0.9

frq <- omega_norm(omeg, Storativ/Transmiss, z = dep, invert = TRUE)
resp <- well_response(frq, Transmiss, Storativ, 
                      Vw, Rs, Ku, B)

test_that("class is correct",{
  expect_is(resp,'wrsp')
})

test_that("method dispatch appropriate",{
  
  expect_is(as.data.frame(resp),'data.frame')
  
  expect_output(expect_is(print(resp),'wrsp'))
  expect_output(expect_message(print(resp)))
  
  expect_is(summary(resp),'summary.wrsp')
  expect_output(expect_is(print(summary(resp)),'summary.wrsp'))
  expect_output(expect_message(print(summary(resp))))
  
})