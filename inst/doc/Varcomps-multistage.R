## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PracTools)

## -----------------------------------------------------------------------------
 require(PracTools)
 data(MDarea.popA)
 BW2stageSRS(MDarea.popA$y1, psuID=MDarea.popA$PSU)
 BW2stageSRS(MDarea.popA$y1, psuID=MDarea.popA$SSU)

## -----------------------------------------------------------------------------
   trtBG <- 10*MDarea.popA$TRACT + MDarea.popA$BLKGROUP
   BW2stageSRS(MDarea.popA$y1, psuID=MDarea.popA$TRACT)
   BW2stageSRS(MDarea.popA$y1, psuID=trtBG)

## -----------------------------------------------------------------------------
 pp.PSU <- table(MDarea.popA$PSU) / nrow(MDarea.popA)
 pp.SSU <- table(MDarea.popA$SSU) / nrow(MDarea.popA)
 BW2stagePPS(MDarea.popA$y1, pp=pp.PSU, psuID=MDarea.popA$PSU)
 BW2stagePPS(MDarea.popA$y1, pp=pp.SSU, psuID=MDarea.popA$SSU)

## -----------------------------------------------------------------------------
 pp.trt <- table(MDarea.popA$TRACT) / nrow(MDarea.popA)
 pp.BG <- table(trtBG) / nrow(MDarea.popA)
 BW2stagePPS(MDarea.popA$y1, pp=pp.trt, psuID=MDarea.popA$TRACT)
 BW2stagePPS(MDarea.popA$y1, pp=pp.BG, psuID=trtBG)

## -----------------------------------------------------------------------------
 M <- length(unique(MDarea.popA$TRACT))
 trtBG <- 10*MDarea.popA$TRACT + MDarea.popA$BLKGROUP
 pp.trt <- rep(1/M,M)
 BW3stagePPS(X=MDarea.popA$y1, pp=pp.trt,
       psuID=MDarea.popA$TRACT, ssuID=trtBG)

## -----------------------------------------------------------------------------
 trtBG <- 10*MDarea.popA$TRACT + MDarea.popA$BLKGROUP
 pp.trt <- table(MDarea.popA$TRACT) / nrow(MDarea.popA)
 BW3stagePPS(X=MDarea.popA$y1, pp=pp.trt,
       psuID=MDarea.popA$TRACT, ssuID=trtBG)

