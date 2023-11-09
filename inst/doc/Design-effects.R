## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
## libraries needed for example 
library(PracTools) 
library(sampling) 
library(ggplot2) 

## Use PracTools smho.N874 data set 
data(smho.N874) 
## Remove hosp.type == 4 as it is out-patient and has no BEDS 
smho <- smho.N874[smho.N874$hosp.type != 4, ] 
## Use co-variate BEDS for MOS 
## Re-code BEDS to have a minimum MOS of 5 
smho$BEDS[smho$BEDS <= 5] <- 5 
## Create 1-draw probability vector based on BEDS 
smho$pi1 <- inclusionprobabilities(smho$BEDS, 1) 
## Create vector for sampling n=50 
pik <- inclusionprobabilities(smho$BEDS, 50) 
## Create sample 
seed <- 20230802 
set.seed(seed) 
sample    <- UPrandomsystematic(pik = pik) 
smho.samp <- smho[sample == 1, ] 
## Create vector of weights 
wgt       <- 1/pik[sample == 1] 

