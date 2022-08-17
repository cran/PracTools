nAuditMUS <- function(MUSVar,        ## Variable used for Monetary Unit Sampling
                      Value.sw,      ## Switch for Positive, Negative, or Absolute Values only
                      CL = 0.90,     ## Confidence Level between 0 and 1; Default is 90%
                      Error.sw,      ## Absolute or Percentage
                      Tol.Error,     ## Tolerable Error Amount/Rate
                      Exp.Error = 0  ## Expected Error Amount/Rate
                         ){

  ## Confirm MUSVar is numeric
  if(is.numeric(MUSVar) == FALSE)
    stop("MUS variable must be numeric.\n")

  ## Confirm MUSVar has no missing values
  if(any(is.na(MUSVar)) == TRUE)
    stop("MUS variable has missing values, which are not allowed.\n")

  ## Give warning where MUSVar has zero-dollar values
  if(any(MUSVar == 0) == TRUE)
    warning("Records with $0 value exist.\n")

  ## Confirm Value Switch is "Pos", "Neg", or "Abs"
  Value.sw <- substr(Value.sw, 1, 3)
  if(Value.sw != "Pos" & Value.sw != "Neg" & Value.sw != "Abs")
    stop("First three characters of Value.sw must be Pos, Neg, or Abs\n")

  ## Confirm Confidence Level is between 0 and 1
  if(CL <= 0 | CL >= 1)
    stop("Confidence level must be between 0 and 1.\n")

  ## Confirm Error Switch is "Amount" or "Percent" (can be abbreviated "Amt" and "Pct")
  if(Error.sw != "Amount"  & Error.sw != "Amt" &
     Error.sw != "Percent" & Error.sw != "Pct")
    stop("Error.sw must be Amount, Amt, Percent, or Pct\n")

  ## Confirm Tolerable Error and Expected Error are between 0 and sum(MUSVar)
  ## when Error.sw is Amount
  if(Error.sw %in% c("Amount", "Amt") & Tol.Error <= 0)
    stop("Tolerable Error must be greater than $0.\n")
  if(Error.sw %in% c("Amount", "Amt") & Exp.Error <  0)
    stop("Expected Error must be greater or equal to $0.\n")
  if(Error.sw %in% c("Amount", "Amt") & Tol.Error >= sum(MUSVar))
    stop("Tolerable Error must be less than the total of MUSvar.\n")
  if(Error.sw %in% c("Amount", "Amt") & Exp.Error >= sum(MUSVar))
    stop("Expected Error must be less than the total of MUSvar.\n")

  ## Confirm Tolerable Error and Expected Error are between 0 and 100
  ## when Error.sw is Percent
  if(Error.sw %in% c("Percent", "Pct") & (Tol.Error < 0 | Tol.Error >= 100))
    stop("Tolerable Error Rate must be in the range [0, 100).\n")
  if(Error.sw %in% c("Percent", "Pct") & (Exp.Error < 0 | Exp.Error >= 100))
    stop("Expected Error Rate must be in the range [0, 100).\n")

  ## Check whether Exp.Error <= Tol.Error
  if(Error.sw %in% c("Amount", "Amt") & Exp.Error >= Tol.Error)
    stop("Exp.Error must be less than Tol.Error")

  ## Create warning message for Error.sw is Pct or Percent, and
  ## Tol.Error or Exp.Error is between 0 and 1 (e.g. .20)
  if(Error.sw %in% c("Percent", "Pct") & (Tol.Error >= 0 & Tol.Error <= 1))
    warning("The errors are in percents. Did you mean to enter a Tol.Error of <= 1%?\n")
  if(Error.sw %in% c("Percent", "Pct") & (Exp.Error >  0 & Exp.Error <  1))
    warning("The errors are in percents. Did you mean to enter a Exp.Error of <= 1%?\n")


  ## Create Population Count and Dollar values
  N <- if(Value.sw == "Pos") sum(MUSVar > 0) else
       if(Value.sw == "Neg") sum(MUSVar < 0) else
                             sum(MUSVar > 0) + sum(MUSVar < 0)

  D <- if(Value.sw == "Pos") sum(MUSVar[MUSVar > 0]) else
       if(Value.sw == "Neg") sum(abs(MUSVar)[MUSVar < 0]) else
                             sum(abs(MUSVar))

  ## Convert Tolerable Error and Expected Error to proportions
  Tol.Error <- if(Error.sw  %in% c("Percent", "Pct")) Tol.Error/100 else
                  Tol.Error/D
  Exp.Error <- if(Error.sw  %in% c("Percent", "Pct")) Exp.Error/100 else
                  Exp.Error/D


  ## MUS Sample size function
  if(Exp.Error == 0) n <- nAuditAttr(TolRate = Tol.Error,
                                     AccDev  = 0,
                                     CL      = CL,
                                     N       = N)$Sample.Size.Hypergeometric

  ## Method for Exp.Error > 0:
  ## Step 1: Identify where hypergeometric sample size / acceptable deviations is greater
  ## than the expected error rate
  ## Step 2: Do linear interpolation as this is most conservative. The formula is:
  ## n = n0 + (Exp.Error - AccDev0/n0)/(AccDev1/n1 - AccDev0/n0)*(n1 - n0)

  if(Exp.Error > 0) n <- {
    if(N <= 100) AccDev <- 0:N else AccDev <- 0:100
    n.sam  <- vector(length = length(AccDev))
    for(i in 1:length(AccDev)){
      n.sam[i]  <- nAuditAttr(TolRate = Tol.Error, CL = CL, N = N,
                              AccDev  = AccDev[i])$Sample.Size.Hypergeometric
    }
    Sample.Exp.Error.Ratio <- AccDev/n.sam
    x0 <- max(Sample.Exp.Error.Ratio[Sample.Exp.Error.Ratio <  Exp.Error], na.rm = TRUE)
    x1 <- min(Sample.Exp.Error.Ratio[Sample.Exp.Error.Ratio >= Exp.Error], na.rm = TRUE)
    n0 <- max(n.sam[Sample.Exp.Error.Ratio <  Exp.Error],                  na.rm = TRUE)
    n1 <- min(n.sam[Sample.Exp.Error.Ratio >= Exp.Error],                  na.rm = TRUE)
    n  <- ceiling(n0 + (Exp.Error - x0)/(x1 - x0)*(n1 - n0))

  }

  structure(list(
       "Value.Range"          = Value.sw,
       "Error.Type"           = Error.sw,
       "Tol.Error.Rate"       = 100*Tol.Error,
       "Exp.Error.Rate"       = 100*Exp.Error,
       "Number.Records"       = N,
       "Sample.Size"          = n,
       "Number.HighVal"       = sum(abs(MUSVar) > D/n),
       "Positive.Pop.Dollars" = D,
       "Conf.Level"           = CL,
       "Sampling.Interval"    = D/n), class = "power.htest")
}
