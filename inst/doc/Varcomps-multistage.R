## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PracTools)

## -----------------------------------------------------------------------------
 require(PracTools)
 data(MDarea.pop)
 BW2stageSRS(MDarea.pop$y1, psuID=MDarea.pop$PSU)
 BW2stageSRS(MDarea.pop$y1, psuID=MDarea.pop$SSU)

## -----------------------------------------------------------------------------
   trtBG <- 10*MDarea.pop$TRACT + MDarea.pop$BLKGROUP
   BW2stageSRS(MDarea.pop$y1, psuID=MDarea.pop$TRACT)
   BW2stageSRS(MDarea.pop$y1, psuID=trtBG)

## -----------------------------------------------------------------------------
 pp.PSU <- table(MDarea.pop$PSU) / nrow(MDarea.pop)
 pp.SSU <- table(MDarea.pop$SSU) / nrow(MDarea.pop)
 BW2stagePPS(MDarea.pop$y1, pp=pp.PSU, psuID=MDarea.pop$PSU)
 BW2stagePPS(MDarea.pop$y1, pp=pp.SSU, psuID=MDarea.pop$SSU)

## -----------------------------------------------------------------------------
 pp.trt <- table(MDarea.pop$TRACT) / nrow(MDarea.pop)
 pp.BG <- table(trtBG) / nrow(MDarea.pop)
 BW2stagePPS(MDarea.pop$y1, pp=pp.trt, psuID=MDarea.pop$TRACT)
 BW2stagePPS(MDarea.pop$y1, pp=pp.BG, psuID=trtBG)

## -----------------------------------------------------------------------------
 M <- length(unique(MDarea.pop$TRACT))
 trtBG <- 10*MDarea.pop$TRACT + MDarea.pop$BLKGROUP
 pp.trt <- rep(1/M,M)
 BW3stagePPS(X=MDarea.pop$y1, pp=pp.trt,
       psuID=MDarea.pop$TRACT, ssuID=trtBG)

## -----------------------------------------------------------------------------
 trtBG <- 10*MDarea.pop$TRACT + MDarea.pop$BLKGROUP
 pp.trt <- table(MDarea.pop$TRACT) / nrow(MDarea.pop)
 BW3stagePPS(X=MDarea.pop$y1, pp=pp.trt,
       psuID=MDarea.pop$TRACT, ssuID=trtBG)

