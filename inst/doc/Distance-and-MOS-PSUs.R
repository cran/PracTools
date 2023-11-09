## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, results='hide'----------------------------------------------------

library(PracTools)
data(Test_Data_US)
GeoDistPSU(Test_Data_US$lat, Test_Data_US$long, "miles", 100, Input.ID = Test_Data_US$ID)


## ----plot.and.histogram, fig.height = 5, fig.width = 7, fig.align = "left"----

g <- GeoDistPSU(Test_Data_US$lat, Test_Data_US$long, "miles", 100, Input.ID = Test_Data_US$ID)

plot(g$PSU.Info$PSU.Mean.Longitude, 
     g$PSU.Info$PSU.Mean.Latitude,
     pch  = 19,
     main = "Plot of PSU Centers",
     xlab = "Longitude",
     ylab = "Latitude")
grid(col = "grey40")

hist(g$PSU.Info$PSU.Max.Dist,
     main = "Histogram of Maximum Within-PSU Distance",
     xlab = "Distance",
     ylab = "Frequency")


## ----Add.PSU------------------------------------------------------------------

## Add PSU from GeoDistPSU
Test_Data_US$psuID <- g$PSU.ID$psuID
## Update GeoDistPSUs for a threshold measure of size of 0.80
m <- GeoDistMOS(lat = Test_Data_US$lat, long = Test_Data_US$long, psuID = Test_Data_US$psuID, n = 15, MOS.var = Test_Data_US$Amount, MOS.takeall = 0.80, Input.ID = Test_Data_US$ID)


## ----MOS.histogram, fig.height = 5, fig.width = 7, fig.align = "left"---------

hist(m$PSU.Max.MOS.Info$psuID.prob,
       breaks = seq(0, 1, 0.05), 
       main = "Histogram of PSU Inclusion Probabilities (Certainties = 1)",
       xlab = "Inclusion Probability",
       ylab = "Frequency")


## ----US.map, warning=FALSE, fig.height = 5, fig.width = 7, fig.align = "left"----
#library(sp)
#library(usmap)
#library(ggplot2)
## Transform PSUs into usmap projection
#g.map  <- cbind(long = g$PSU.Info$PSU.Mean.Longitude, 
#          lat  = g$PSU.Info$PSU.Mean.Latitude) 
#g.map  <- as.data.frame(g.map)
#g.proj <- usmap::usmap_transform(g.map,
#             input_names  = c("long", "lat"),
#             output_names = c("Long", "Lat"))
#plot_usmap(color = "gray") +
#        geom_point(data = g.proj, aes(x = Long, y = Lat))

## ----BW2stageSRS--------------------------------------------------------------

BW2stageSRS(X = Test_Data_US$Y, psuID = Test_Data_US$psuID, lonely.SSU = "zero")


## ----BW2stagePPS--------------------------------------------------------------

pp <- tapply(Test_Data_US$Amount, Test_Data_US$psuID, sum)/sum(Test_Data_US$Amount)
BW2stagePPS(X = Test_Data_US$Y, pp = pp, psuID = Test_Data_US$psuID, lonely.SSU = "zero")


## ----MOS.PSU.merge, message=FALSE---------------------------------------------

library(dplyr)
Test_Data_US <- Test_Data_US %>% mutate(ID = as.character(ID))
Test_Data_US <- inner_join(Test_Data_US, m$PSU.ID.Max.MOS, by=c("ID" = "Input.ID"))


## ----BW2stageSRS.MOS----------------------------------------------------------

BW2stageSRS(X = Test_Data_US$Y, psuID = Test_Data_US$psuID.new, lonely.SSU = "zero")


## ----Certainties.MOS----------------------------------------------------------

certs <- (1:nrow(m$PSU.Max.MOS.Info))[m$PSU.Max.MOS.Info$psuID.prob > 0.8]
certID <- m$PSU.Max.MOS.Info[certs, "psuID.new"]
certID


## ----Onedraw.MOS--------------------------------------------------------------
## Create vector of 1-draw probabilities
pp <- tapply(Test_Data_US$Amount, Test_Data_US$psuID.new, sum)/sum(Test_Data_US$Amount)
## Subset Test_Data_US for non-certainties
sub.Test_Data_US <- Test_Data_US[!(Test_Data_US$psuID.new %in% certID),]
## Subset vector of 1-draw probabilities for non-certainties
sub.pp <- pp[-certs]

## ----BW2stagePPS.MOS----------------------------------------------------------
## Rescale sub.pp to sum to 1
sub.pp <- sub.pp/sum(sub.pp)
BW2stagePPS(X = sub.Test_Data_US$Y, pp = sub.pp, psuID = sub.Test_Data_US$psuID.new, 
            lonely.SSU = "zero")


