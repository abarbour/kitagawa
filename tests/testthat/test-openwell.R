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

test_that("wang == hsieh when leak = 0",{
  
  Transmiss <- 1e-6
  Storativ <- 1e-6
  omeg <- 10**seq(-3,3, by=0.1) # Hz
  dep <- 30 # meters
  frq <- omega_norm(omeg, Storativ/Transmiss, z = dep, invert = TRUE)
  
  
  hsieh <- open_well_response(frq, Transmiss, Storativ, z=dep, model = 'hsieh')
  wang  <- open_well_response(frq, Transmiss, Storativ, z=dep, model = 'wang', leak = 0)

  expect_equal(hsieh$Response, wang$Response)
  
})


test_that("wang spot check phase",{

  
  Transmiss <- 1e-6
  Storativ  <- 1e-2
  omeg      <- 1.9322736 / 86400 # Hz
  leak      <- c(10^(-10), 10^(-3.094118))
  
  wang  <- open_well_response(omeg, Transmiss, Storativ,
                              z=dep, 
                              Rs. = 0.1,
                              model = 'wang', 
                              leak = leak,
                              freq.units = 'Hz', 
                              as.pressure = FALSE)
  
  # Check result based on digitized values
  expect_equal(as.numeric(Arg(wang$Response[,2]) * 180/pi), 
                    c(-46.658323, 78.372966),
                    tolerance = 0.5, scale = 1)
  
  
  Transmiss <- 1e-4
  Storativ  <- 1e-2
  omeg      <- 1.9322736 / 86400 # Hz
  leak      <- c(10^(-10.600000), 10^(-2))
  
  wang  <- open_well_response(omeg, Transmiss, Storativ,
                              z=dep, 
                              Rs. = 0.1,
                              model = 'wang', 
                              leak = leak,
                              freq.units = 'Hz', 
                              as.pressure = FALSE)
  
  # Check result based on digitized value
  expect_equal(as.numeric(Arg(wang$Response[,2]) * 180/pi), 
                    c(-1.6020025, 89.6370463),
               tolerance = 0.3, scale = 1)
  
  
  Transmiss <- 10
  Storativ  <- 1e-2
  omeg      <- 1.9322736 / 86400 # Hz
  leak      <- c(10^(-20), 10^(0))
  
  wang  <- open_well_response(omeg, Transmiss, Storativ,
                              z=dep, 
                              Rs. = 0.1,
                              model = 'wang', 
                              leak = leak,
                              freq.units = 'Hz', 
                              as.pressure = FALSE)
  
  # Check result end members
  expect_equivalent(as.numeric(Arg(wang$Response[,2]) * 180/pi), 
                    c(0, 90), tolerance = 0.1, scale = 1)
  
  
})


