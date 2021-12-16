## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
#  devtools::install_github("LandSciTech/LSTDConnect",
#                           ref = "paper")

## ----load_libs----------------------------------------------------------------
library(LSTDConnect)
library(raster)
library(samc)
library(microbenchmark)
library(ggplot2)

## ----get_data-----------------------------------------------------------------
res <- LSTDConnect::ghm
plot(res)

## ----get_occ------------------------------------------------------------------
occ <- 101 - res
plot(occ)

## ----get_mort-----------------------------------------------------------------
high_mort <- 0.25
low_mort <- 0.01
cutoff <- 80

mort <- res
mort[mort >= cutoff] <- high_mort
mort[mort < cutoff] <- low_mort

## ----step_1-------------------------------------------------------------------
samc_cache <- LSTDConnect::samc(resistance = res, absorption = mort, 
                                directions = 4)

## ----step_2-------------------------------------------------------------------
dists <- LSTDConnect::distribution(samc = samc_cache, occ = occ, time = 1000)
plot(dists$occ)

## ----samc_samc----------------------------------------------------------------
tr_list <- list(fun = function(x) 1 / mean(x),
                dir = 8, 
                sym = TRUE)

samc_cache_comp <- samc::samc(data = res, absorption = mort, tr_args = tr_list)
dists_comp <- samc::distribution(samc = samc_cache_comp, occ = occ, time = 1000)
dists_mapped <- map(samc_cache_comp, dists_comp)
plot(dists_mapped)

## ----comp---------------------------------------------------------------------
diff <- dists_mapped - dists$occ
plot(diff)

## ----comp_hist----------------------------------------------------------------
hist(values(diff), main = "Error")

## -----------------------------------------------------------------------------
t <- 1000

samc_obj <- samc::samc(data = res,
                       absorption = mort,
                       tr_args = tr_list)  
samc_obj_custom <- LSTDConnect::samc(resistance = res, 
                                     absorption = mort, 
                                     directions = 8)

mbm <- microbenchmark(
  "samc" = {
    short_disp <- samc::distribution(samc = samc_obj,
                                     occ = mort,
                                     time = t)
    gc()
  },
  "LSTDConnect" = {
    short_disp_custom <- LSTDConnect::distribution(samc = samc_obj_custom,
                                                   occ = occ, 
                                                   time = t)
    gc()
  },
  times = 100,
  unit = "s") 
mbm

## ----benchmark_plot-----------------------------------------------------------
plot(mbm)

