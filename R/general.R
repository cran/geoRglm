"cut0.calc.mixed.gauss" <- 
  function(value, parms)
{
  if(ncol(parms$var)==1) parms$var <- as.vector(parms$var)
  if(any(parms$var == 0)) parms$var[parms$var == 0] <- 1e-12
  temp <- array(pnorm((value - parms$mean)/sqrt(parms$var)), dim = c(nrow(parms$mean), ncol(parms$mean)))
  temp2 <- apply(temp, 1, mean)
  return(temp2)
}

"cut0.calc.mixed.t" <- 
  function(value, parms, df)
{
  if(ncol(parms$var)==1) parms$var <- as.vector(parms$var)
  if(any(parms$var == 0)) parms$var[parms$var == 0] <- 1e-12
  temp <- array(pt((value - parms$mean)/sqrt((parms$var * (df - 2))/df), df = df), dim = c(nrow(parms$mean), ncol(parms$mean)))
  temp2 <- apply(temp, 1, mean)
  return(temp2)
}

"model.glm.control" <- 
  function(trend.d = "cte", trend.l = "cte", cov.model = "matern", kappa = 0.5, aniso.pars = NULL, lambda = 0)
{
  cov.model <-
    match.arg(cov.model,
              choices = c("matern", "exponential", "gaussian",
                "spherical", "circular", "cubic", "wave", "power",
                "powered.exponential", "cauchy", "gneiting",
                "gneiting.matern", "pure.nugget"))
  if(cov.model == "powered.exponential" & (kappa <= 0 | kappa > 2))
    stop("krige.bayes: for power exponential correlation model the parameter kappa must be in the interval \(0,2\]")
  ##  if(any(cov.model == c("exponential", "gaussian", "spherical",
  ##           "circular", "cubic", "wave", "powered.exponential",
  ##           "cauchy", "gneiting", "pure.nugget")))
  ##    kappa <- NULL
  if(!is.null(aniso.pars)){ 
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be a vector with two elements: rotation angle (in radians) and anisotropy ratio (a number > 1)")
  }
  res <- list(trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, kappa = kappa, aniso.pars = aniso.pars, lambda = lambda)
  class(res) <- "model.geoRglm"
  return(res)
}


"multgauss" <- 
  function(cov)
{
  if(is.R())
    return(crossprod(chol(cov), rnorm(n=ncol(cov))))
  else
    return(rmvnorm(ncol(cov), cov = cov))
}

"phi.number" <- 
  function(phi.prior = c("uniform", "exponential", "fixed", "squared.reciprocal", "reciprocal"))
{
  ## WARNING: codes above must be the same as in the C function logprior_phi
  ## 
  phi.prior <- match.arg(phi.prior)
  phinumber <- switch(phi.prior,
                      uniform = as.integer(1),
                      exponential = as.integer(2),
                      fixed = as.integer(3),
                      squared.reciprocal = as.integer(4),
                      reciprocal = as.integer(5),
                      stop("wrong or no specification of phi.prior"))
  return(phinumber)
}


"mcmc.control" <- 
  function(S.scale, Htrunc="default", S.start, burn.in=0, thin=10, n.iter=1000*thin, phi.start="default",  phi.scale=NULL)
{
  if(missing(S.scale)) stop("S.scale parameter must to be provided for MCMC-proposal")
  if(missing(S.start)) S.start<-"default"
  if(is.numeric(S.start))
    if(length(S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
  else   
    if(!any(S.start == "default")) stop(" S.start must be a vector of same dimension as data ")              
  res <- list(S.scale = S.scale, Htrunc = Htrunc, S.start = S.start, burn.in = burn.in, thin = thin, n.iter = n.iter,
              phi.start = phi.start, phi.scale=phi.scale)
  class(res) <- "mcmc.geoRglm"
  return(res)
}


"prior.glm.control" <- 
  function(beta.prior = c("flat", "normal", "fixed"), beta = NULL, beta.var.std = NULL, 
           sigmasq.prior = c("uniform", "sc.inv.chisq", "reciprocal", "fixed"), sigmasq = NULL, df.sigmasq = NULL,
           phi.prior = c("uniform","exponential", "fixed", "squared.reciprocal", "reciprocal"), phi = NULL, 
           phi.discrete = NULL, tausq.rel = 0)
{
  beta.prior <- match.arg(beta.prior)
  sigmasq.prior <- match.arg(sigmasq.prior)
  phi.prior <- match.arg(phi.prior)
  if(beta.prior == "fixed" & is.null(beta))
    stop("argument \"beta\" must be provided with fixed value for this parameter")
  if(beta.prior == "normal"){
    if(is.null(beta) | is.null(beta.var.std))
      stop("arguments \"beta\" and \"beta.var.std\" must be provided when using normal prior for the parameter beta")
    if((length(beta))^2 != length(beta.var.std))
      stop(" beta and beta.var.std have incompatible dimensions")
    if(any(beta.var.std != t(beta.var.std)))
      stop(" non symmetric matrix in beta.var.std")
    if(inherits(try(chol(beta.var.std)), "try-error"))
          stop(" matrix in beta.var.std is not positive definit")
  }
  ##
  if(sigmasq.prior == "fixed" & is.null(sigmasq))
    stop("argument \"sigmasq\" must be provided when the parameter sigmaq is fixed")
  if(sigmasq.prior == "sc.inv.chisq")
    if(is.null(sigmasq) | is.null(df.sigmasq))
      stop("arguments \"sigmasq\" and \"df.sigmasq\" must be provided for inverse chisq prior")
  if(!is.null(sigmasq))
    if(sigmasq < 0) stop("negative values not allowed for \"sigmasq\"")
  if(sigmasq.prior == "reciprocal"){
    warning("This choice of sigmasq.prior gives an improper posterior !!!!!!! \n")
    sigmasq <- 0
    df.sigmasq <- 0
  }
  if(sigmasq.prior == "uniform"){
    sigmasq <- 0
    df.sigmasq <- -2
  }
  ##
  if(phi.prior == "fixed"){
    if(is.null(phi)){
      stop("argument \"phi\" must be provided with fixed prior for this parameter")
    }
  }
  else{
    if(phi.prior == "exponential" & (is.null(phi) | (length(phi) > 1)))
      stop("argument \"phi\" must be provided when using the exponential prior for the parameter phi")
    if(any(phi.prior == c("reciprocal", "squared.reciprocal")) & any(phi.discrete == 0)){
      warning("degenerated prior at phi = 0. Excluding value phi.discrete[1] = 0")
      phi.discrete <- phi.discrete[phi.discrete > 1e-12]
    }
    if(!is.null(phi.discrete)){
      discrete.diff <- diff(phi.discrete)
      if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
        stop("The current implementation requires equally spaced values in the argument \"phi.discrete\"\n")
    }
    if(phi.prior != "exponential") phi <- NULL
  }
  if(any(phi.discrete < 0))
    stop("negative values not allowed for parameter phi")
  ##
  if(is.null(tausq.rel)) stop("argument \"tausq.rel\" must be provided")
  ##
  ip <- list(beta=list(), sigmasq=list(), phi=list())
  ##
  if(beta.prior == "fixed"){
    ip$beta$status <- "fixed"
    ip$beta$fixed.value <- beta 
  }
  else{
    ip$beta <- list(dist = beta.prior)
    if(beta.prior == "flat")
      ip$beta$pars <- c(0, +Inf)
    if(beta.prior == "normal"){
      if(length(beta) == 1)
        ip$beta$pars <- c(mean=beta, var.std=beta.var.std)
      else
        ip$beta$pars <- list(mean=beta, var.std=beta.var.std)
    }
  }
  ##
  if(sigmasq.prior == "fixed"){
    ip$sigmasq$status <- "fixed"
    ip$sigmasq$fixed.value <- sigmasq 
  }
  else{
    ip$sigmasq <- list(dist = sigmasq.prior)
    if(sigmasq.prior == "reciprocal")
      ip$sigmasq$pars <- c(df=0, var=+Inf)
    if(sigmasq.prior == "uniform")
      ip$sigmasq$pars <- c(df=-2, var=+Inf)
    if(sigmasq.prior == "sc.inv.chisq")
      ip$sigmasq$pars <- c(df=df.sigmasq, var=sigmasq)
  }
  ##
  if(phi.prior == "fixed"){
    ip$phi$status <- "fixed"
    ip$phi$fixed.value <- phi
  }
  else{
    ip$phi$dist <- phi.prior
    if(is.null(phi.discrete))
      stop("phi.discrete must be given when parameter phi is random")
    else{
      ip$phi$probs <- switch(phi.prior,
                             uniform = rep(1/length(phi.discrete), length(phi.discrete)),
                             exponential = (1/phi) * exp(- phi.discrete/phi),
                             squared.reciprocal = (1/(phi.discrete^2))/sum(1/(phi.discrete^2)),
                             reciprocal = (1/phi.discrete)/sum(1/phi.discrete))
      names(ip$phi$probs) <- phi.discrete
    }
    if(phi.prior == "exponential") ip$phi$pars <- c(ip$phi$pars, exp.par=phi)
  }
  ##
  ip$tausq.rel <- list(status = "fixed", fixed.value = tausq.rel)
  ##
  res <- list(beta.prior = beta.prior, beta = beta, beta.var.std = beta.var.std, sigmasq.prior = sigmasq.prior, sigmasq = sigmasq, 
              df.sigmasq = df.sigmasq, phi.prior = phi.prior, phi = phi, 
              phi.discrete = phi.discrete, tausq.rel = tausq.rel, priors.info = ip)
  class(res) <- "prior.geoRglm"
  return(res)
}

"BCinv.aux" <- 
  function(z,lambda)
{
  return(BCtransform(z,lambda = lambda, inverse=TRUE)$data)
}

"output.glm.control" <-
  function(sim.posterior, sim.predict, keep.mcmc.sim, quantile, threshold, inference, messages.screen)
{
  ##
  ## Assigning default values
  ##
  if(missing(sim.posterior)) sim.posterior <- TRUE
  if(missing(sim.predict)) sim.predict <- FALSE
  if(missing(keep.mcmc.sim)) keep.mcmc.sim <- TRUE
  if(missing(quantile))  quantile.estimator <- NULL
  else quantile.estimator <- quantile
  if(missing(threshold)) probability.estimator <- NULL 
  else probability.estimator <- threshold
  if(missing(inference)) inference <- TRUE
  if(missing(messages.screen)) messages.screen <- TRUE
  ##
  ##
  if(!is.null(quantile.estimator)){
    if(is.numeric(quantile.estimator))
      if(any(quantile.estimator) < 0 | any(quantile.estimator) > 1)
        stop("quantiles indicators must be numbers in the interval [0,1]\n")
    if(any(quantile.estimator == TRUE)) quantile.estimator <- c(0.025, 0.5, 0.975)
    if(!inference) {
      warning("prediction not performed; quantile.estimator is set to NULL \n")
      quantile.estimator <- NULL
    }
  }
  if(!is.null(probability.estimator)){
    if(!is.numeric(probability.estimator))
      stop("threshold must be a numeric value (or vector) of cut-off value(s)\n")
    if(!inference) {
      warning("prediction not performed; probabilitites above a threshold cannot be computed\n")
      probability.estimator <- NULL
    }
  }
  res <- list(sim.posterior = sim.posterior, sim.predict = sim.predict, keep.mcmc.sim = keep.mcmc.sim,
              quantile.estimator = quantile.estimator, probability.estimator = probability.estimator,
              inference = inference, messages.screen = messages.screen)
  class(res) <- "output.geoRglm"
  return(res)
}
