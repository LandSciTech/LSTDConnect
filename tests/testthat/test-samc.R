# Prepare

library(raster)
res <- LSTDConnect::ghm
occ <- 101 - res
high_mort <- 0.25
low_mort <- 0.01
cutoff <- 80
mort <- res
mort[mort >= cutoff] <- high_mort
mort[mort < cutoff] <- low_mort

test_that("samc works", {
  samc_cache <- LSTDConnect::samc(resistance = res, absorption = mort,
                                  directions = 8)
  dists <- LSTDConnect::distribution(samc = samc_cache, occ = occ, time = 10)
  dists <- raster::as.matrix(dists$occ)
  
  tr_list <- list(fun = function(x) 1 / mean(x),
                  dir = 8, 
                  sym = TRUE)
  
  samc_cache_comp <- suppressWarnings(
    samc::samc(data = res, absorption = mort, tr_args = tr_list))
  dists_comp <- samc::distribution(samc = samc_cache_comp, occ = occ, time = 10)
  dists_mapped <- samc::map(samc_cache_comp, dists_comp)
  dists_mapped <- raster::as.matrix(dists_mapped)
  
  expect_equal(dists, dists_mapped, tolerance = 1e-15)
})
