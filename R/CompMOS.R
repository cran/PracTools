CompMOS <- function(dsn           = NULL,  ## Data frame used for Composite MOS calculations
                    psuID         = NULL,  ## PSU Cluster ID
                    n.PSU         = NULL,  ## PSU sample size
                    domain        = NULL,  ## Vector of domain variable names
                    domain.req.n  = NULL,  ## Required sample from each domain
                    exp.domain.rr = NULL   ## Expected response rate for each domain 
                                           ## in percentage between 0 and 1
                    ) {

  ## Confirm dsn is data frame
  if(is.data.frame(dsn) == FALSE)
    stop("dsn must be a data frame.\n")
  
  ## Confirm psuID is not null
  if(is.vector(psuID) == FALSE)
    stop("psuID must be a vector.\n")
  
  ## Confirm domain is not null
  if(is.vector(domain) == FALSE)
    stop("domain must be a vector of the domain variables.\n")
  
  ## Confirm exp.domain.rr is not null
  if(is.vector(exp.domain.rr) == FALSE)
    stop("exp.domain.rr must be a vector of the expected domain response rates.\n")
  
  ## Confirm domain.req.n is not null
  if(is.vector(domain.req.n) == FALSE)
    stop("domain.req.n must be a vector of the required sample size in each domain.\n")
  
  ## Confirm sample size is positive number
  if(is.null(n.PSU) == TRUE)
    stop("n.PSU (PSU Sample Size) must be provided.\n")
  if(is.numeric(n.PSU) == FALSE)
    stop("n.PSU (PSU Sample Size) must be numeric.\n")
  if(n.PSU <= 0)
    stop("n.PSU (PSU Sample Size) must be greater than 0\n")

  ## Calculate Domain populations
  Nh <- colSums(dsn[, domain])

  ## Calculate Domain expected sample size
  nh <- domain.req.n / exp.domain.rr 
  
  ## Calculate Domain sampling fraction
  fh <- nh/Nh
  
  ## Calculate Total Workload and PSU Workload
  wrkld     <- sum(nh)
  PSU.wrkld <- sum(nh)/n.PSU
    
  ## Calculate Composite Measure of Size
  ## Create data frame replicating fh for each psuID
  m <- do.call("rbind", replicate(n = nrow(dsn), fh, simplify = FALSE))
  ## Create data frame for Composite MOS calculation
  ## Note: trick is to have psuID and domain read in separately as cbind
  CompMOS.dsn <- cbind(psuID, dsn[, c(domain)])
  ## Calculate Composite MOS for each psuID
  CompMOS.dsn$CompMOS = rowSums(m * CompMOS.dsn[, -1])
  
  ## Calculate PSU probability of inclusion
  CompMOS.dsn$probi <- n.PSU*CompMOS.dsn$CompMOS/sum(CompMOS.dsn$CompMOS)
    
  ## Calculate within-PSU Domain probabilities and sample sizes, and perform 
  ## a feasibility check on pi(k|i)
  ## Create data frame 
  CompMOS.domain    <- CompMOS.dsn[, domain]
  ## Create PSU within-domain sampling fractions
  CompMOS.domain.sf             <- PSU.wrkld * m / CompMOS.dsn$CompMOS
  colnames(CompMOS.domain.sf)   <- paste0(colnames(CompMOS.domain.sf),".samp.frac")
  ## Create PSU within-domain sample sizes
  CompMOS.domain.ss             <- CompMOS.domain * PSU.wrkld * m / CompMOS.dsn$CompMOS
  colnames(CompMOS.domain.ss)   <- paste0(colnames(CompMOS.domain.ss),".samp.size")
  ## Perform PSU within-domain sampling feasibility checks
  CompMOS.domain.piki           <- (m <= CompMOS.dsn$CompMOS/PSU.wrkld)
  colnames(CompMOS.domain.piki) <- paste0(colnames(CompMOS.domain.piki),".feas.chk")
    
  warn <- "None"
  if(any(CompMOS.domain.piki==FALSE)){
    warning("One or more PSUs fails the feasibility check for domain sampling. \n")
    warn <- "One or more PSUs fails the feasibility check for domain sampling."
  }
    
  ## Create output data frames  
  ## Create PSU-level output data frame
  CompMOS.psuID <- cbind(CompMOS.dsn, 
                         CompMOS.domain.sf, 
                         CompMOS.domain.ss,
                         CompMOS.domain.piki)
    
  ## Create Composite MOS information data frame
  CompMOS.design           <- t(rbind(Nh, exp.domain.rr, domain.req.n, nh, fh))
  colnames(CompMOS.design) <- c("Dom.Pop", "Dom.E(RR)", "Req.Dom.Sample", 
                                "Req.Adj.Dom.Sample", "Dom.Samp.Fraction")
  
  ## Create Composite MOS operations data frame
  CompMOS.Ops <- cbind(n.PSU, wrkld, PSU.wrkld)
  colnames(CompMOS.Ops) <- c("Number.PSUs", "Sample.Workload", "PSU.Workload")
 
  ## Output Composite MOS data for PSU, Design, and Ops
  out <- list(warning        = warn,  
              CompMOS.psuID  = CompMOS.psuID,
              CompMOS.design = CompMOS.design,
              CompMOS.Ops    = CompMOS.Ops
              )
  
  ## Return output
  return(out)
}
