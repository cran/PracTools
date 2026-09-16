pclass <- function (formula, data=NULL, link="logit", numcl=5, type, design=NULL, seed=NULL, rng.preds = 0.05){
    if (!(link %in% c("logit", "probit", "cloglog")))
        stop("link must be logit, probit, or cloglog.\n")
    if (!(type %in% c("wtd", "unwtd")))
        stop("type must be wtd or unwtd.\n")
  
  
    if (type == "unwtd"){
        if (is.null(data)) stop("With type=unwtd data cannot be NULL.\n")
        reg <- glm(formula, family = binomial(link = link), data)
    }
    else {
        if (is.null(design))
            stop("design must be specified for weighted analysis.\n")
        if (link == "logit")
            reg <- survey::svyglm(formula, design, family = quasibinomial(link = "logit"))
        if (link == "probit")
            reg <- survey::svyglm(formula, design, family = quasibinomial(link = "probit"))
        if (link == "cloglog")
            reg <- survey::svyglm(formula, design, family = quasibinomial(link = "cloglog"))
    }
    L.hat <- reg$linear.predictors

    if (link == "logit"){
        preds <- exp(L.hat) / (1 + exp(L.hat) )
    }
    if (link == "probit"){
        preds <- pnorm(L.hat)
    }
    if (link == "cloglog"){
        preds <- 1- exp(-exp(L.hat) )
    }
    
    if (max(preds) - min(preds) <= rng.preds){
      p.class <- rep(1, length(preds))
      preds <- rep(mean(preds), length(preds))
      warning("Range of response propensities is less than or equal to rng.preds. A single adjustment class is formed.\n")
    }
    else{
      if (!is.null(seed)) {
        if (!is.null(seed)) {
          if (exists(".Random.seed", envir = .GlobalEnv)) {
            old.seed <- get(".Random.seed", envir = .GlobalEnv)
            on.exit(assign(".Random.seed", old.seed, envir = .GlobalEnv), add = TRUE)
          }
          set.seed(seed)
        }
      }
      preds <- jitter(preds)
      pcuts <- seq(0, 1, 1/numcl)
      quintiles <- quantile(preds, probs = pcuts)
      p.class <- cut(preds, breaks = quintiles, include.lowest=TRUE)
    }
    
    list(p.class = p.class,
         propensities = preds)
}
