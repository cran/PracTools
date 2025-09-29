## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Ex1, echo=TRUE, eval=TRUE------------------------------------------------

library(PracTools)
library(kableExtra)
data(hospital)
HOSP <- PracTools::hospital
HOSP$sqrt.x <- sqrt(HOSP$x)
sam <- sample(nrow(HOSP), 50)
N1.resp <- HOSP[sam, ]
N2.nonresp  <- HOSP[-sam, ]

## Create lm object using "known" data; no intercept model
lm.obj  <- lm(y ~ 0 + sqrt.x + x, data = N1.resp)

## Create range of values to use as delta for difference in means
delta <- mean(HOSP$y) - mean(HOSP$y) * seq(.6, 1, by=0.05)

## Run SampStop function and output to object S
S <- SampStop(lm.obj  = lm.obj,
              formula = ~ 0 + sqrt.x + x,
              n1.data = N1.resp, 
              yvar    = "y", 
              n2.data = N2.nonresp, 
              p       = seq(0.2, 0.6, by=0.05), 
              delta   = delta, 
              seed    = .Random.seed[413]) 

kableExtra::kable(S$Input,  caption = "SampStop Input")
kableExtra::kable(head(S$Output, n=15),
      caption = "SampStop Output: First 15 Observations")

## ----Ex1.plot, echo=TRUE, eval=TRUE, fig.align="center", fig.width = 7--------

library(ggplot2)
## Convert S to data frame
S1 <- as.data.frame(S$Output)

## Create factor category over probability of response and number of responders
p.nresp <- paste(S1$`Pr(response)`, S1$`Exp no. resps`, sep=", ")

ggplot(S1, aes(x = `diff in means`, 
              y = `Pr(smaller diff)`, 
              colour = factor(p.nresp))) +
  geom_point() +
  geom_line(linewidth=1.1) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) + 
  labs(title = "Probability of Response by Delta",
       x = "delta", 
       y = "Pr(|e1 - e2|<= delta)", 
       colour = "Pr(Resp), Number \nof Responders") 


