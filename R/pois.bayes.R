
"mcmc.bayes.pois.log" <- 
  function(data, units.m, trend, mcmc.input, cov.model, kappa, tausq.rel, distcoords, ss.sigma, df, phi.prior, phi.discrete, phi)
{
  ##
#### This is the MCMC engine for the Bayesian analysis of a spatial Poisson log Normal model
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5 
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n) Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, log(data), 0) - log(units.m))
  }
  else S <- as.vector(mcmc.input$S.start)
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(phi.prior == "fixed"){
    phi.discrete <- phi
    phi.prior <- "uniform"
  }
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1) {
      p.discr.rejects <- 1 - 2 * pnorm((phi.discrete[nmphi] - phi.discrete[1])/(2 * (nmphi - 1)), sd = sqrt(phi.scale), lower.tail = FALSE )
      if(p.discr.rejects > 0.2)
         warning(paste("probability of rejecting phi-proposal because of discr. is \n", p.discr.rejects))
    }
  }
  ##                                                                      
##### ---------- sampling ----------- ###### 
  cov.model.number <- cor.number(cov.model)
  phi.prior.number <- phi.number(phi.prior)
  if(is.null(phi)) {
    if(phi.prior == "exponential")
      stop("Must give the parameter in exponential prior-distribution for phi")
    else e.mean <- 0
  }
  else e.mean <- phi
   if(is.vector(trend)) beta.size <- 1
  else beta.size <- ncol(trend)
  n.sim <- floor(n.iter/thin)
  Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  phi.sample <- as.double(rep(phi, n.sim))
  if(is.R()){
    result <-  .C("mcmcrun4",
                  as.integer(n),
                  as.double(data),
                  as.double(units.m),
                  as.double(as.vector(t(trend))),
                  as.integer(beta.size),
                  as.integer(cov.model.number),
                  as.double(kappa),
                  as.double(tausq.rel),
                  as.double(distcoords),
                  as.double(S.scale),
                  as.double(phi.scale),
                  as.double(Htrunc),
                  as.integer(n.iter),
                  as.integer(thin),
                  as.integer(burn.in),
                  as.double(ss.sigma),
                  as.integer(df),
                  as.integer(phi.prior.number),
                  as.double(phi.discrete),
                  as.integer(nmphi),
                  as.double(e.mean),
                  Sdata = Sdata,
                  phi.sample = phi.sample, DUP=FALSE)[c("Sdata", "phi.sample")]
  }
  else{
    result <- .C("mcmcrun4",
                 COPY = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                 as.integer(n),
                 as.double(data),
                 as.double(units.m),
                 as.double(as.vector(t(trend))),
                 as.integer(beta.size),
                 as.integer(cov.model.number),
                 as.double(kappa),
                 as.double(tausq.rel),
                 as.double(distcoords),
                 as.double(S.scale),
                 as.double(phi.scale),
                 as.double(Htrunc),
                 as.integer(n.iter),
                 as.integer(thin),
                 as.integer(burn.in),
                 as.double(ss.sigma),
                 as.integer(df),
                 as.integer(phi.prior.number),
                 as.double(phi.discrete),
                 as.integer(nmphi),
                 as.double(e.mean),
                 Sdata = Sdata,
                 phi.sample = phi.sample)[c("Sdata", "phi.sample")]
  }
  attr(result$Sdata, "dim") <- c(n, n.sim)
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  return(result)
}


"mcmc.bayes.pois.boxcox" <- 
  function(data, units.m, trend, mcmc.input, cov.model, kappa, tausq.rel, distcoords, ss.sigma, df, phi.prior, phi.discrete, phi, lambda)  
{
  ##
#### This is the MCMC engine for the Bayesian analysis of a spatial Poisson boxcox Normal model, when lambda >0
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, ((data/units.m)^lambda -1)/lambda, 0) )         
  }
  else S <- as.vector(mcmc.input$S.start)
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5 
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n) Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(phi.prior == "fixed"){
    phi.discrete <- phi
    phi.prior <- "uniform"
  }
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1) {
      p.discr.rejects <- 1 - 2 * pnorm((phi.discrete[nmphi] - phi.discrete[1])/(2 * (nmphi - 1)), sd = sqrt(phi.scale),  lower.tail = FALSE)
      if(p.discr.rejects > 0.2)
         warning(paste("probability of rejecting phi-proposal because of discr. is \n", p.discr.rejects))
    }
  }
  ##                                                                      
##### ---------- sampling ----------- ###### 
  cov.model.number <- cor.number(cov.model)
  phi.prior.number <- phi.number(phi.prior)
  if(is.null(phi)) {
    if(phi.prior == "exponential")
      stop("Must give the parameter in exponential prior-distribution for phi")
    else e.mean <- 0
  }
  else e.mean <- phi
  if(is.vector(trend)) beta.size <- 1
  else beta.size <- ncol(trend)
  n.sim <- floor(n.iter/thin)
  Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  phi.sample <- as.double(rep(phi, n.sim))
  if(is.R()){
    result <-  .C("mcmcrun4boxcox",
                  as.integer(n),
                  as.double(data),
                  as.double(units.m),              
                  as.double(as.vector(t(trend))),
                  as.integer(beta.size),
                  as.integer(cov.model.number),
                  as.double(kappa),
                  as.double(tausq.rel),
                  as.double(distcoords),
                  as.double(S.scale),
                  as.double(phi.scale),
                  as.double(Htrunc),
                  as.integer(n.iter),
                  as.integer(thin),
                  as.integer(burn.in),
                  as.double(ss.sigma),
                  as.integer(df),
                  as.integer(phi.prior.number),
                  as.double(phi.discrete),
                  as.integer(nmphi),
                  as.double(e.mean),
                  as.double(lambda),
                  Sdata = Sdata,
                  phi.sample = phi.sample, DUP=FALSE)[c("Sdata", "phi.sample")]    
  }
  else{
    result <- .C("mcmcrun4boxcox",
                 COPY = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                 as.integer(n),
                 as.double(data),
                 as.double(units.m),                 
                 as.double(as.vector(t(trend))),
                 as.integer(beta.size),
                 as.integer(cov.model.number),
                 as.double(kappa),
                 as.double(tausq.rel),
                 as.double(distcoords),
                 as.double(S.scale),
                 as.double(phi.scale),
                 as.double(Htrunc),
                 as.integer(n.iter),
                 as.integer(thin),
                 as.integer(burn.in),
                 as.double(ss.sigma),
                 as.integer(df),
                 as.integer(phi.prior.number),
                 as.double(phi.discrete),
                 as.integer(nmphi),
                 as.double(e.mean),
                 as.double(lambda),
                 Sdata = Sdata,
                 phi.sample = phi.sample)[c("Sdata", "phi.sample")]    
  }
  attr(result$Sdata, "dim") <- c(n, n.sim)
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  return(result)
}


"mcmc.bayes.conj.pois.log" <- 
  function(data, units.m, meanS, ttvbetatt, mcmc.input, cov.model, kappa, tausq.rel, distcoords, ss.sigma, df, phi.prior,
           phi.discrete, phi)
{
  ##
  ## This is the MCMC engine for the Bayesian analysis (with normal prior for beta) of a spatial binomial logit Normal model
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$S.start == "default")) {
    S <- as.vector(ifelse(data > 0, log(data), 0) - log(units.m))
  }
  else S <- as.vector(mcmc.input$S.start)
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5 
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n) Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(phi.prior == "fixed"){
    phi.discrete <- phi
    phi.prior <- "uniform"
  }
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1) {
      p.discr.rejects <- 1 - 2 * pnorm((phi.discrete[nmphi] - phi.discrete[1])/(2 * (nmphi - 1)), sd = sqrt(phi.scale),  lower.tail = FALSE)
      if(p.discr.rejects > 0.2)
        warning(paste("probability of rejecting phi-proposal because of discr. is \n", p.discr.rejects))
    }
  }
  ##                                                                      
  ## ---------- sampling ----------- ###### 
  cov.model.number <- cor.number(cov.model)
  phi.prior.number <- phi.number(phi.prior)
  if(is.null(phi)) {
    if(phi.prior == "exponential")
      stop("Must give the parameter in exponential prior-distribution for phi")
    else e.mean <- 0
  }
  else e.mean <- phi
  if(is.null(ttvbetatt)) ttvbetatt <- matrix(0,beta.size,beta.size)
  n.sim <- floor(n.iter/thin)
  Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  phi.sample <- as.double(rep(phi, n.sim)) 
  if(is.R()){
    result <-  .C("mcmcrun5",
                  as.integer(n),
                  as.double(data),
                  as.double(units.m),
                  as.double(as.vector(meanS)),
                  as.double(as.vector(ttvbetatt)),
                  as.integer(cov.model.number),
                  as.double(kappa),
                  as.double(tausq.rel),
                  as.double(distcoords),
                  as.double(S.scale),
                  as.double(phi.scale),
                  as.double(Htrunc),
                  as.integer(n.iter),
                  as.integer(thin),
                  as.integer(burn.in),
                  as.double(ss.sigma),
                  as.integer(df),
                  as.integer(phi.prior.number),
                  as.double(phi.discrete),
                  as.integer(nmphi),
                  as.double(e.mean),
                  Sdata = Sdata,
                  phi.sample = phi.sample, DUP=FALSE)[c("Sdata", "phi.sample")]
  }
  else{
    result <- .C("mcmcrun5",
                 COPY = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                 as.integer(n),
                 as.double(data),
                 as.double(units.m),
                 as.double(as.vector(meanS)),
                 as.double(as.vector(ttvbetatt)),
                 as.integer(cov.model.number),
                 as.double(kappa),
                 as.double(tausq.rel),
                 as.double(distcoords),
                 as.double(S.scale),
                 as.double(phi.scale),
                 as.double(Htrunc),
                 as.integer(n.iter),
                 as.integer(thin),
                 as.integer(burn.in),
                 as.double(ss.sigma),
                 as.integer(df),
                 as.integer(phi.prior.number),
                 as.double(phi.discrete),
                 as.integer(nmphi),
                 as.double(e.mean),
                 Sdata = Sdata,
                 phi.sample = phi.sample)[c("Sdata", "phi.sample")]
  }
  attr(result$Sdata, "dim") <- c(n, n.sim)
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  return(result)
}

"mcmc.bayes.conj.pois.boxcox" <- 
  function(data, units.m, meanS, ttvbetatt, mcmc.input, cov.model, kappa, tausq.rel, distcoords, ss.sigma, df, phi.prior,
           phi.discrete, phi, lambda)
{
  ##
  ## This is the MCMC engine for the Bayesian analysis (with normal prior for beta) of a spatial binomial logit Normal model
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$S.start == "default")) {
    S <- as.vector(ifelse(data > 0, ((data/units.m)^lambda -1)/lambda, 0) )  
  }
  else S <- as.vector(mcmc.input$S.start)
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5 
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n) Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  if(phi.prior == "fixed"){
    phi.discrete <- phi
    phi.prior <- "uniform"
  }
  if(any(mcmc.input$phi.start=="default")) phi <- median(phi.discrete)
  else  phi <- mcmc.input$phi.start
  nmphi <-  length(phi.discrete)
  if(is.null(mcmc.input$phi.scale)) {
    if(nmphi > 1) stop("mcmc.input$phi.scale not given ")
    else phi.scale <- 0
  }
  else {
    phi.scale <- mcmc.input$phi.scale
    if(nmphi > 1) {
      p.discr.rejects <- 1 - 2 * pnorm((phi.discrete[nmphi] - phi.discrete[1])/(2 * (nmphi - 1)), sd = sqrt(phi.scale), lower.tail = FALSE )
      if(p.discr.rejects > 0.2)
        warning(paste("probability of rejecting phi-proposal because of discr. is \n", p.discr.rejects))
    }
  }
  ##                                         
  ## ---------- sampling ----------- ###### 
  cov.model.number <- cor.number(cov.model)
  phi.prior.number <- phi.number(phi.prior)
  if(is.null(phi)) {
    if(phi.prior == "exponential")
      stop("Must give the parameter in exponential prior-distribution for phi")
    else e.mean <- 0
  }
  else e.mean <- phi
  if(is.null(ttvbetatt)) ttvbetatt <- matrix(0,beta.size,beta.size)
  n.sim <- floor(n.iter/thin)
  Sdata <- as.double(as.vector(c(S, rep(0, (n.sim - 1) * n))))
  phi.sample <- as.double(rep(phi, n.sim)) 
  if(is.R()){
    result <-  .C("mcmcrun5boxcox",
                  as.integer(n),
                  as.double(data),
                  as.double(units.m),
                  as.double(as.vector(meanS)),
                  as.double(as.vector(ttvbetatt)),
                  as.integer(cov.model.number),
                  as.double(kappa),
                  as.double(tausq.rel),
                  as.double(distcoords),
                  as.double(S.scale),
                  as.double(phi.scale),
                  as.double(Htrunc),
                  as.integer(n.iter),
                  as.integer(thin),
                  as.integer(burn.in),
                  as.double(ss.sigma),
                  as.integer(df),
                  as.integer(phi.prior.number),
                  as.double(phi.discrete),
                  as.integer(nmphi),
                  as.double(e.mean),
                  as.double(lambda),
                  Sdata = Sdata,
                  phi.sample = phi.sample, DUP=FALSE)[c("Sdata", "phi.sample")]
  }
  else{
    result <- .C("mcmcrun5boxcox",
                 COPY = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
                   FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
                 as.integer(n),
                 as.double(data),
                 as.double(units.m),
                 as.double(as.vector(meanS)),
                 as.double(as.vector(ttvbetatt)),
                 as.integer(cov.model.number),
                 as.double(kappa),
                 as.double(tausq.rel),
                 as.double(distcoords),
                 as.double(S.scale),
                 as.double(phi.scale),
                 as.double(Htrunc),
                 as.integer(n.iter),
                 as.integer(thin),
                 as.integer(burn.in),
                 as.double(ss.sigma),
                 as.integer(df),
                 as.integer(phi.prior.number),
                 as.double(phi.discrete),
                 as.integer(nmphi),
                 as.double(e.mean),
                 as.double(lambda),
                 Sdata = Sdata,
                 phi.sample = phi.sample)[c("Sdata", "phi.sample")]
  }
  attr(result$Sdata, "dim") <- c(n, n.sim)
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  return(result)
}


"pois.krige.bayes" <- 
  function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", locations = "no", 
           model, prior, mcmc.input, output)
{
###########
  if(missing(geodata))
    geodata <- list(coords=coords, data=data)
  if(is.R()) require(mva)
  call.fc <- match.call()
  seed <- .Random.seed
  do.prediction <- ifelse(all(locations == "no"), FALSE, TRUE)
  ##
  ## Checking data configuration
  ##
  if(is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  if(nrow(coords) != length(data))
    stop("number of data is different of number of data locations (coordinates)")
  n <- length(data)
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n)
  }
  ####
  ## reading model input
  ##
  if(missing(model))
    model <- model.glm.control()
  else{
    if(is.null(class(model)) || class(model) != "model.geoRglm"){
      if(!is.list(model))
        stop("pois.krige.bayes: the argument model only takes a list or an output of the function model.glm.control")
      else{
        model.names <- c("trend.d", "trend.l", "cov.model", "kappa", "aniso.pars", "lambda") 
        model.user <- model
        model <- list()
        if(length(model.user) > 0){
          for(i in 1:length(model.user)){
            n.match <- match.arg(names(model.user)[i], model.names)
            model[[n.match]] <- model.user[[i]]
          }
        }    
        if(is.null(model$trend.d)) model$trend.d <- "cte"  
        if(is.null(model$trend.l)) model$trend.l <- "cte"  
        if(is.null(model$cov.model)) model$cov.model <- "matern"  
        if(is.null(model$kappa)) model$kappa <- 0.5
        if(is.null(model$lambda)) model$lambda <- 0
        model <- model.glm.control(trend.d = model$trend.d,
                                   trend.l = model$trend.l,
                                   cov.model = model$cov.model,
                                   kappa = model$kappa,
                                   aniso.pars = model$aniso.pars,
                                   lambda = model$lambda)
      }
    }
  }
  cov.model <- model$cov.model
  kappa <- model$kappa
  tausq.rel <- prior$tausq.rel
  lambda <- model$lambda
  if(lambda < 0) stop ("lambda < 0 is not allowed")
  ## reading prior input
  ##
  if(missing(prior)) stop("pois.krige.bayes: argument prior must be given")
  else{
    if(is.null(class(prior)) || class(prior) != "prior.geoRglm"){
      if(!is.list(prior))
        stop("pois.krige.bayes: the argument prior only takes a list or an output of the function prior.glm.control")
      else{
        prior.names <- c("beta.prior", "beta", "beta.var.std", "sigmasq.prior",
                         "sigmasq", "df.sigmasq", "phi.prior", "phi", "phi.discrete", "tausq.rel") 
        prior.user <- prior
        prior <- list()
        if(length(prior.user) > 0){
          for(i in 1:length(prior.user)){
            n.match <- match.arg(names(prior.user)[i], prior.names)
            prior[[n.match]] <- prior.user[[i]]
          }
        } 
        if(is.null(prior$beta.prior)) prior$beta.prior <- "flat"
        if(is.null(prior$sigmasq.prior)) prior$sigmasq.prior <- "uniform"
        if(is.null(prior$phi.prior)) prior$phi.prior <- "uniform"
        if(is.null(prior$tausq.rel)) prior$tausq.rel <- 0
        prior <- prior.glm.control(beta.prior = prior$beta.prior,
                                   beta = prior$beta, beta.var.std = prior$beta.var.std,
                                   sigmasq.prior = prior$sigmasq.prior,
                                   sigmasq = prior$sigmasq,  df.sigmasq = prior$df.sigmasq,
                                   phi.prior = prior$phi.prior,
                                   phi = prior$phi, phi.discrete = prior$phi.discrete, 
                                   tausq.rel = prior$tausq.rel)
      }
    }
  }
  beta.prior <- prior$beta.prior
  beta <- prior$beta
  beta.var <- prior$beta.var.std
  sigmasq.prior <- prior$sigmasq.prior
  if(sigmasq.prior == "fixed") sigmasq <- prior$sigmasq
  else{
    df.sigmasq <- prior$df.sigmasq
    S2.prior <- prior$sigmasq
  }
  phi.prior <- prior$phi.prior 
  phi <- prior$phi
  if(sigmasq.prior == "fixed" & phi.prior != "fixed") stop("option for fixed sigmasq and random phi not implemented")
  if(phi.prior != "fixed") phi.discrete <- prior$phi.discrete
  else phi.discrete <- phi
  ##
  ## reading output options
  ##
  if(missing(output))
    output <- output.glm.control()
  else{
    if(is.null(class(output)) || class(output) != "output.geoRglm"){
      if(!is.list(output))
        stop("pois.krige.bayes: the argument output only takes a list or an output of the function output.glm.control")
      else{
        output.names <- c("sim.posterior","sim.predict", "keep.mcmc.sim","quantile","threshold","inference","messages.screen")      
        output.user <- output
        output <- list()
        if(length(output.user) > 0){
          for(i in 1:length(output.user)){
            n.match <- match.arg(names(output.user)[i], output.names)
            output[[n.match]] <- output.user[[i]]
          }
        }
        if(is.null(output$sim.posterior)) output$sim.posterior <- TRUE
        if(is.null(output$sim.predict)) output$sim.predict <- FALSE
        if(is.null(output$keep.mcmc.sim)) output$keep.mcmc.sim <- TRUE 
        if(is.null(output$inference)) output$inference <- TRUE
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE        
        output <- output.glm.control(sim.posterior = output$sim.posterior,
                                     sim.predict = output$sim.predict,
                                     keep.mcmc.sim = output$keep.mcmc.sim, quantile = output$quantile,
                                     threshold = output$threshold, inference = output$inference,
                                     messages.screen = output$messages.screen)
      }
    }
  }
  quantile.estimator <- output$quantile.estimator
  probability.estimator <- output$probability.estimator  
  inference <- output$inference
  messages.screen <- output$messages.screen
  ## check == here
  data.dist <- as.vector(dist(coords))
  if(round(1000000000000. * min(data.dist)) == 0) stop("Two coords are identical; not allowed.")
  ##
  trend.d <- model$trend.d
  if(messages.screen) {
    if(inherits(trend.d, "formula")) cat("pois.krige.bayes: analysis using covariates provided by the user\n")
    else {
      if(trend.d == "cte") cat("pois.krige.bayes: analysis assuming constant mean\n")
      if(trend.d == "1st") cat("pois.krige.bayes: analysis assuming a 1st degree polynomial trend\n")
      if(trend.d == "2nd") cat("pois.krige.bayes: analysis assuming a 2nd degree polynomial trend\n")
    }
  }
  trend.data <- trend.spatial(trend = trend.d, geodata = geodata)
  dimnames(coords) <- list(NULL, NULL)
  dimnames(trend.data) <- list(NULL, NULL)
  beta.size <- ncol(trend.data)
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  if(beta.prior == "normal" |  beta.prior == "fixed"){
    if(beta.size != length(beta))
      stop("pois.krige.bayes: size of beta incompatible with the trend model (covariates)")
  }
  aniso.pars <- model$aniso.par
  if(!is.null(aniso.pars)) coords.transf <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
  else coords.transf <- coords
  ##
  ## checking prediction locations
  ##
  if((inference) & (do.prediction)){
    ## Checking the consistency between coords, locations, and trends
    trend.l <- model$trend.l
    if(is.vector(locations)) {
      if(length(locations) == 2) {
        locations <- t(as.matrix(locations))
        warning("only one location to be predicted (in two-dimensional space) \n")
      }
      else locations <- as.matrix(cbind(locations, 0))
    }
    else locations <- as.matrix(locations)
    ni <- nrow(locations)
    if(is.null(trend.l)) stop("trend.l needed for prediction")
    if(inherits(trend.d, "formula") | inherits(trend.l, "formula")) {
      if(!(inherits(trend.d, "formula")) | !(inherits(trend.l, "formula")))
        stop("trend.d and trend.l must have similar specification")
    }
    else if(trend.d != trend.l) stop("trend.l is different from trend.d")
    kb.results <- list(posterior = list(), predictive = list())
  }
  else {
    if(do.prediction & messages.screen) cat(paste("need to specify inference=TRUE to make predictions \n"))
    kb.results <- list(posterior = list(), predictive = paste("prediction not performed"))
    do.prediction <- FALSE
  }
  ##
  ## ##### preparing for MCMC -------------------------------------------------------
  ##
  if(missing(mcmc.input)) stop("pois.krige.bayes: argument mcmc.input must be given")
  else{
    if(is.null(class(mcmc.input)) || class(mcmc.input) != "mcmc.geoRglm"){
      if(!is.list(mcmc.input))
        stop("pois.krige.bayes: the argument mcmc.input only takes a list or an output of the function mcmc.control")
      else{
        mcmc.input.names <- c("S.scale", "Htrunc", "S.start", "burn.in", "thin", "n.iter", "phi.start",  "phi.scale")    
        mcmc.input.user <- mcmc.input
        mcmc.input <- list()
        if(length(mcmc.input.user) > 0){
          for(i in 1:length(mcmc.input.user)){
            n.match <- match.arg(names(mcmc.input.user)[i], mcmc.input.names)
            mcmc.input[[n.match]] <- mcmc.input.user[[i]]
          }
        }
        if(is.null(mcmc.input$Htrunc)) mcmc.input$Htrunc <- "default"
        if(is.null(mcmc.input$S.start)) mcmc.input$S.start <- "default"
        if(is.null(mcmc.input$burn.in)) mcmc.input$burn.in <- 0
        if(is.null(mcmc.input$thin)) mcmc.input$thin <- 10
        if(is.null(mcmc.input$n.iter)) mcmc.input$n.iter <- 1000*mcmc.input$thin
        if(is.null(mcmc.input$phi.start)) mcmc.input$phi.start <- "default"
        mcmc <- mcmc.control(S.scale = mcmc.input$S.scale,Htrunc=mcmc.input$Htrunc,S.start=mcmc.input$S.start,
                             burn.in=mcmc.input$burn.in,thin=mcmc.input$thin,n.iter=mcmc.input$n.iter,
                             phi.start=mcmc.input$phi.start,phi.scale=mcmc.input$phi.scale)
      }
    }
  }
  ##
  if(beta.prior == "fixed" | beta.prior == "normal") mean.d <- as.vector(trend.data %*% beta)
  else mean.d <- 0
  if(sigmasq.prior != "fixed"){
    if(beta.prior == "flat") df.model <- n - beta.size + df.sigmasq
    else df.model <- n + df.sigmasq
  }
  else df.model <- Inf
  if(beta.prior == "normal"){
    if(beta.size > 1) ttvbetatt <- trend.data%*%beta.var%*%t(trend.data)
    else ttvbetatt <- trend.data%*%t(trend.data)*beta.var
  }  
  else ttvbetatt <- NULL
  if(sigmasq.prior == "fixed") {     ### implies that phi is fixed !
    if(beta.prior == "normal"){
      QQt <- chol(varcov.spatial(coords = coords.transf, cov.model = cov.model, kappa = kappa, cov.pars = c(1,phi),
                                 nugget = tausq.rel)$varcov + ttvbetatt )
    }
    else {
      QQt <- varcov.spatial(coords = coords.transf, cov.model = cov.model, kappa = kappa, cov.pars = c(1,phi),
                            nugget = tausq.rel, only.decomposition = TRUE, try.another.decomposition = FALSE)$sqrt.varcov
    }
    if(beta.prior == "flat"){
      if(is.R()) QQtinvtt <- crossprod(solve(QQt),trend.data)
      else QQtinvtt <- crossprod(solve.upper(QQt),trend.data)
      Dmat <- diag(n) - QQtinvtt %*% solve.geoR(crossprod(QQtinvtt)) %*% t(QQtinvtt)
    }
    else Dmat <- NULL
  }
  ##
############----------PART 2 ------------##############################
############-----------MCMC -------------##############################
  ##
  if(sigmasq.prior == "fixed"){ 
    if(lambda == 0){
      gauss.post <- mcmc.pois.log(data = data, units.m = units.m, meanS = mean.d, QQ = sqrt(sigmasq)*t(QQt), mcmc.input = mcmc.input,
                                Dmat = Dmat)
    }
    else{
      gauss.post <- mcmc.pois.boxcox(data = data, units.m = units.m, meanS = mean.d, QQ = sqrt(sigmasq)*t(QQt), mcmc.input = mcmc.input,
                                   lambda=lambda,  Dmat = Dmat)
    }
  } 
  else {
    kb.results$posterior$phi <- list()
    ## take care re-using gauss.post !
    if(beta.prior == "flat"){
      if(lambda == 0){ 
        gauss.post <- mcmc.bayes.pois.log(data=data, units.m=units.m, trend=trend.data, mcmc.input=mcmc.input, cov.model=cov.model, 
                                        kappa=kappa, tausq.rel = tausq.rel, distcoords=distdiag(coords.transf), 
                                        ss.sigma = df.sigmasq*S2.prior, df = df.model, phi.prior = phi.prior,
                                        phi.discrete = phi.discrete, phi = phi)
      }
      else{
        gauss.post <- mcmc.bayes.pois.boxcox(data=data, units.m=units.m, trend=trend.data, mcmc.input=mcmc.input, cov.model=cov.model, 
                                           kappa=kappa, tausq.rel = tausq.rel, distcoords=distdiag(coords.transf),
                                           ss.sigma = df.sigmasq*S2.prior, df = df.model,  phi.prior = phi.prior,
                                             phi.discrete = phi.discrete, phi = phi, lambda = lambda)
      }
    }
    else{
      if(lambda == 0){ 
        gauss.post <- mcmc.bayes.conj.pois.log(data=data, units.m=units.m, meanS = mean.d, ttvbetatt = ttvbetatt, mcmc.input=mcmc.input,
                                               cov.model=cov.model,  kappa=kappa, tausq.rel = tausq.rel,
                                               distcoords=distdiag(coords.transf),  ss.sigma = df.sigmasq*S2.prior, df = df.model,
                                               phi.prior = phi.prior, phi.discrete = phi.discrete, phi = phi)
      }
      else{ 
        gauss.post <- mcmc.bayes.conj.pois.boxcox(data=data, units.m=units.m, meanS = mean.d, ttvbetatt = ttvbetatt,
                                                  mcmc.input=mcmc.input, cov.model=cov.model,  kappa=kappa, tausq.rel = tausq.rel,
                                                  distcoords=distdiag(coords.transf),  ss.sigma = df.sigmasq*S2.prior, df = df.model,
                                                  phi.prior = phi.prior, phi.discrete = phi.discrete, phi = phi, lambda = lambda)
      }       
    }
    kb.results$posterior$phi$sample <- gauss.post$phi.sample
    gauss.post <- gauss.post$Sdata
  }  
  ## ######### removing junk #########################################
  if(is.R()){
    if(sigmasq.prior == "fixed"){
      remove(list = c("QQt", "Dmat"))
      if(beta.prior == "flat") remove("QQtinvtt")
    }
    remove("ttvbetatt") 
  }
  else {
    if(sigmasq.prior == "fixed"){
      remove(list = c("QQt", "Dmat"), frame = sys.nframe())
      if(beta.prior == "flat") remove("QQtinvtt", frame = sys.nframe())
    }
    remove("ttvbetatt", frame = sys.nframe())
  }  
  ##           
##############-------------PART 3----------######################
##############------------prediction-------######################
  ## 
  n.sim <- ncol(gauss.post)
  if(inference) {
    temp.post <- list()
    temp.post$beta.mean <- array(NA, dim = c(beta.size, n.sim))
    temp.post$beta.var <- array(NA, dim = c(beta.size, beta.size, n.sim))
    temp.post$S2 <- rep(0, n.sim)
    if(do.prediction) {
      temp.pred <- list()
      temp.pred$mean <- array(NA, dim = c(ni, n.sim))
      temp.pred$var <- array(NA, dim = c(ni, n.sim))
      if(output$sim.predict) {
        num.pred <- 1
        kb.results$predictive$simulations <- array(NA, dim = c(ni, n.sim))
      }
      else {
        num.pred <- 0
        kb.results$predictive$simulations <- " no simulations from the preditive distribution "
      }
    }
    else num.pred <- 0
    model.temp <- model
    model.temp$lambda <- 1
    output.temp <- list(n.posterior = 0, n.predictive = num.pred, messages.screen = FALSE)
    prior.temp <- prior
    prior.temp$phi.prior <- "fixed"
    prior.temp$phi.discrete <- NULL
    prior.temp$priors.info <- NULL
    if(phi.prior == "fixed" || length(phi.discrete) == 1) {
      if(phi.prior == "fixed") prior.temp$phi <- phi
      else prior.temp$phi <- phi.discrete
      temp.result <- krige.bayes.extnd(data = gauss.post, coords = coords, locations = locations, 
                                       model = model.temp, prior = prior.temp, output = output.temp)
      temp.post$beta.mean <- temp.result$posterior$beta$pars$mean        
      temp.post$beta.var <- temp.result$posterior$beta$pars$var
      temp.post$S2 <- temp.result$posterior$sigmasq$pars$S2
      if(do.prediction) {
        temp.pred$mean <- temp.result$predictive$mean
        temp.pred$var <- temp.result$predictive$variance
        if(output$sim.predict){
          if(lambda == 0) kb.results$predictive$simulations <- exp(temp.result$predictive$simulations)
          else kb.results$predictive$simulations <- BC.inv(temp.result$predictive$simulations, lambda=lambda)
        }
      }      
    }
    else {    
      len.phi.discrete <- length(phi.discrete)      
      step.phi.discrete <- phi.discrete[2] - phi.discrete[1]
      phi.table <- rep(0,len.phi.discrete)
      for(i in 1:len.phi.discrete){
        phi.table[i] <- sum(ifelse(abs(kb.results$posterior$phi$sample-phi.discrete[i])<0.5*step.phi.discrete,1,0))
      }
      phi.sample.unique <- phi.discrete[phi.table>0]
      phi.table <- phi.table[phi.table>0]
      len.phi.un <- length(phi.sample.unique)
      indic.phi <- array(rep(0, len.phi.un * max(phi.table)), dim = c(len.phi.un, max(phi.table)))
      numbers <- 1:length(kb.results$posterior$phi$sample)
      for(i in 1:len.phi.un) {
        temp.num <- numbers[abs(kb.results$posterior$phi$sample-phi.sample.unique[i])<0.5*step.phi.discrete]
        indic.phi[i, 1:length(temp.num)] <- temp.num
      }
      for(i in 1:len.phi.un){
        if(phi.table[i] == 1){
          prior.temp$phi <- phi.sample.unique[i]
          temp.result <- krige.bayes(data = gauss.post[, indic.phi[i,1]], coords = coords, locations = locations, 
                                     model = model.temp, prior = prior.temp, output = output.temp)
          temp.post$beta.mean[, indic.phi[i,1]] <- temp.result$posterior$beta$pars$mean
          temp.post$beta.var[,  , indic.phi[i,1]] <- temp.result$posterior$beta$pars$var
          temp.post$S2[indic.phi[i,1]] <- temp.result$posterior$sigmasq$pars$S2
          if(do.prediction){        
            temp.pred$mean[, indic.phi[i, 1]] <- temp.result$predictive$mean            
            temp.pred$var[, indic.phi[i, 1]] <- temp.result$predictive$variance            
            if(output$sim.predict) {
              kb.results$predictive$simulations[, indic.phi[i, 1]] <- BC.inv(temp.result$predictive$simulations, lambda)
            }
          }
        }
        else {
          prior.temp$phi <- phi.sample.unique[i]
          temp.result <- krige.bayes.extnd(data = gauss.post[, indic.phi[i,1:phi.table[i]]], coords = coords, locations = locations, 
                                           model = model.temp, prior = prior.temp, output = output.temp)  
          temp.post$beta.mean[, indic.phi[i,1:phi.table[i]]] <- temp.result$posterior$beta$pars$mean         
          temp.post$beta.var[,  , indic.phi[i,1:phi.table[i]]] <- temp.result$posterior$beta$pars$var
          temp.post$S2[indic.phi[i,1:phi.table[i]]] <- temp.result$posterior$sigmasq$pars$S2
          if(do.prediction) {
            temp.pred$mean[, indic.phi[i, 1:phi.table[i]]] <- temp.result$predictive$mean
            temp.pred$var[, indic.phi[i, 1:phi.table[i]]] <- temp.result$predictive$variance          
            if(output$sim.predict) {
              kb.results$predictive$simulations[ , indic.phi[i, 1:phi.table[i]]] <- BC.inv(temp.result$predictive$simulations, lambda)  
            }
          }
        }
      }
    }    
    if(is.R()) remove("temp.result")
    else remove("temp.result", frame = sys.nframe())
    if(do.prediction) {
      ##
      d0mat <- loccoords(coords, locations)
      loc.coincide <- (apply(d0mat < 1e-10, 2, sum) == 1)
      ##
      ##------- calculating medians and uncertainty
      ##
      temp.med <- apply(temp.pred$mean, 1, median)
      temp.unc <- sqrt(apply(temp.pred$mean, 1, var) + apply(temp.pred$var, 1, median))
      not.accurate <- (!loc.coincide)
      diffe <- pmixed(temp.med, temp.pred,df.model)-0.5
      temp.med.new <- temp.med[not.accurate]+0.1*(temp.med[not.accurate]+0.1) # to get started
      inv.sl <- rep(0,ni)
      parms.temp <- list()
      while(any(not.accurate)){
        parms.temp$mean<-temp.pred$mean[not.accurate,,drop=FALSE]
        parms.temp$var<-temp.pred$var[not.accurate,,drop=FALSE]
        diffe.new <- pmixed(temp.med.new, parms.temp,df.model)-0.5
        inv.sl[not.accurate] <- (temp.med.new-temp.med[not.accurate])/(diffe.new-diffe[not.accurate])
        temp.med[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), temp.med.new,temp.med[not.accurate])
        diffe[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), diffe.new, diffe[not.accurate])
        not.accurate[not.accurate] <- ifelse(abs(diffe[not.accurate])>0.0005, TRUE, FALSE)
        temp.med.new <- temp.med[not.accurate] - diffe[not.accurate]*inv.sl[not.accurate]
      }
      temp.upper <- qnorm(rep(0.975, ni), mean = temp.med, sd = temp.unc)
      not.accurate <- (!loc.coincide)
      diffe <- pmixed(temp.upper, temp.pred,df.model)-0.975
      temp.upper.new <- temp.upper[not.accurate]+0.5*(temp.upper[not.accurate]+0.5) # to get started
      inv.sl <- rep(0,ni)      
      while(any(not.accurate)){
        parms.temp$mean<-temp.pred$mean[not.accurate,,drop=FALSE]
        parms.temp$var<-temp.pred$var[not.accurate,,drop=FALSE]
        diffe.new <- pmixed(temp.upper.new, parms.temp,df.model)-0.975
        inv.sl[not.accurate] <- (temp.upper.new-temp.upper[not.accurate])/(diffe.new-diffe[not.accurate])
        temp.upper[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), temp.upper.new,temp.upper[not.accurate])
        diffe[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), diffe.new, diffe[not.accurate])
        not.accurate[not.accurate] <- ifelse(abs(diffe[not.accurate])>0.0005, TRUE, FALSE)
        temp.upper.new <- temp.upper[not.accurate] - diffe[not.accurate]*inv.sl[not.accurate]
      }      
      temp.lower <- qnorm(rep(0.025, ni), mean = temp.med, sd = temp.unc)
      not.accurate <- (!loc.coincide)
      diffe <- pmixed(temp.lower, temp.pred,df.model)-0.025
      temp.lower.new <- temp.lower[not.accurate]+0.5*(temp.lower[not.accurate]+0.5) # to get started
      inv.sl <- rep(0,ni)
      while(any(not.accurate)){
        parms.temp$mean<-temp.pred$mean[not.accurate,,drop=FALSE]
        parms.temp$var<-temp.pred$var[not.accurate,,drop=FALSE]
        diffe.new <- pmixed(temp.lower.new,parms.temp,df.model)-0.025
        inv.sl[not.accurate] <- (temp.lower.new-temp.lower[not.accurate])/(diffe.new-diffe[not.accurate])
        temp.lower[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), temp.lower.new,temp.lower[not.accurate])
        diffe[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), diffe.new, diffe[not.accurate])
        not.accurate[not.accurate] <- ifelse(abs(diffe[not.accurate])>0.0005, TRUE, FALSE)
        temp.lower.new <- temp.lower[not.accurate] - diffe[not.accurate]*inv.sl[not.accurate]
      }
      if(any(loc.coincide)){
        temp.med[loc.coincide] <- apply(temp.pred$mean[loc.coincide,,drop=FALSE], 1, median)
        temp.upper[loc.coincide] <- apply(temp.pred$mean[loc.coincide,,drop=FALSE], 1, quantile, probs = 0.975)
        temp.lower[loc.coincide] <- apply(temp.pred$mean[loc.coincide,,drop=FALSE], 1, quantile, probs = 0.025) 
      }
      kb.results$predictive$median <- BC.inv(temp.med,lambda)
      kb.results$predictive$uncertainty <- (BC.inv(temp.upper,lambda) - BC.inv(temp.lower,lambda))/4
      ## calculating quantiles
      if(is.null(quantile.estimator)) {
        kb.results$predictive$quantiles <- as.data.frame(cbind(BC.inv(temp.lower,lambda), BC.inv(temp.med,lambda), BC.inv(temp.upper,lambda)))
      }
      if(is.numeric(quantile.estimator)) {
        nmq <- length(quantile.estimator)
        if(nmq > 1) {
          temp.quan <- matrix(NA, ni, nmq)
          dig <- rep(3, nmq)
          for(i in 1:nmq) {
            while(quantile.estimator[i] != round(quantile.estimator[i], digits = dig[i])) dig[i] <-dig[i] + 1
            temp.quan[, i] <- qnorm(rep(quantile.estimator[i], ni), mean = temp.med, sd = temp.unc)
            if(any(loc.coincide)) temp.quan[loc.coincide, i] <- temp.med[loc.coincide]
            not.accurate <- (!loc.coincide)
            diffe <- pmixed(temp.quan[,i], temp.pred,df.model)-quantile.estimator[i]
            numb <- 0.1+abs(quantile.estimator[i]-0.5)
            temp.quan.new <- temp.quan[not.accurate,i]+numb*(temp.quan[not.accurate,i]+numb) # to get started
            inv.sl <- rep(0,ni)
            while(any(not.accurate)) {
              parms.temp$mean <-temp.pred$mean[not.accurate,,drop=FALSE]
              parms.temp$var <-temp.pred$var[not.accurate,,drop=FALSE]
              diffe[not.accurate] <- pmixed(temp.quan.new,parms.temp,df.model)-quantile.estimator[i]
              inv.sl[not.accurate] <- (temp.quan.new-temp.quan[not.accurate, i])/(diffe.new-diffe[not.accurate])
              temp.quan[not.accurate, i] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), temp.quan.new,temp.quan[not.accurate, i])
              diffe[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), diffe.new, diffe[not.accurate])
              not.accurate[not.accurate] <- ifelse(abs(diffe[not.accurate])>0.0005, TRUE, FALSE)
              temp.quan.new <- temp.quan[not.accurate, i] - diffe[not.accurate]*inv.sl[not.accurate]
            }            
            if(any(loc.coincide)){
              temp.quan[loc.coincide,i] <- apply(temp.pred$mean[loc.coincide,,drop=FALSE], 1, quantile, probs = quantile.estimator[i])
            }
          }
          kb.results$predictive$quantiles <- as.data.frame(BC.inv(temp.quan,lambda))
        }
        else {
          dig <- 3
          while(quantile.estimator != round(quantile.estimator,digits = dig)) dig <- dig + 1
          temp.quan <- qnorm(rep(quantile.estimator,ni), mean = temp.med, sd = temp.unc)
          not.accurate <- (!loc.coincide)
          diffe <- pmixed(temp.quan, temp.pred,df.model)-quantile.estimator
          numb <- 0.1+abs(quantile.estimator-0.5)
          temp.quan.new <- temp.quan[not.accurate]+numb*(temp.quan[not.accurate]+numb) # to get started
          inv.sl <- rep(0,ni)
          while(any(not.accurate)) {
            parms.temp$mean <-temp.pred$mean[not.accurate,,drop=FALSE]
            parms.temp$var <-temp.pred$var[not.accurate,,drop=FALSE]
            diffe[not.accurate] <- pmixed(temp.quan.new,parms.temp,df.model)-quantile.estimator
            inv.sl[not.accurate] <- (temp.quan.new-temp.quan[not.accurate])/(diffe.new-diffe[not.accurate])
            temp.quan[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), temp.quan.new,temp.quan[not.accurate])
            diffe[not.accurate] <- ifelse(abs(diffe[not.accurate]) > abs(diffe.new), diffe.new, diffe[not.accurate])
            not.accurate[not.accurate] <- ifelse(abs(diffe[not.accurate])>0.0005, TRUE, FALSE)
            temp.quan.new <- temp.quan[not.accurate] - diffe[not.accurate]*inv.sl[not.accurate]
          }
          if(any(loc.coincide)){
            temp.quan[loc.coincide] <- apply(temp.pred$mean[loc.coincide,,drop=FALSE], 1, quantile, probs = quantile.estimator)
          }
          kb.results$predictive$quantiles <- as.vector(BC.inv(temp.quan,lambda))
        }
      }
      if(is.null(quantile.estimator)){
        qname <- rep(0, 3)
        qname[1] <- paste("q0.025", sep = "")
        qname[2] <- paste("q0.5", sep = "")
        qname[3] <- paste("q0.975", sep = "")
        names(kb.results$predictive$quantiles) <- qname
      }
      if(is.numeric(quantile.estimator) && nmq > 1) {
        qname <- rep(0, length(quantile.estimator))
        for(i in 1:length(quantile.estimator))
          qname[i] <- paste("q", 100 * quantile.estimator[i], sep = "")
        names(kb.results$predictive$quantiles) <- qname
      }
      ##
      ## ------ probability estimators
      ##
      if(!is.null(probability.estimator)) {
        if(lambda == 0) transf.probab <- ifelse(probability.estimator > 0, log(probability.estimator), -1e+17)
        else transf.probab <- ifelse(probability.estimator > 0, (probability.estimator^lambda-1)/lambda, -1e+17)
        if(length(transf.probab) == 1) transf.probab <- as.vector(transf.probab)
        for (ii in 1:length(transf.probab)){
          thresh.vec <- rep(transf.probab[ii],ni)
          kb.results$predictive$probability <- round(pmixed(thresh.vec, temp.pred, df.model), digits = 3) 
        }
      }
      if(is.R()) remove("temp.pred")
      else remove("temp.pred", frame = sys.nframe())
      ## 
      if(messages.screen) cat("pois.krige.bayes: Prediction performed \n")
    }
    else {
      kb.results$predictive <- "no locations to perform prediction were provided"
      if(messages.screen) cat(paste("Only Bayesian estimation of model parameters "))
    }
    ##
    ##----- calculating posterior summaries ----------------##
    ##
    if(beta.prior == "fixed") kb.results$posterior$beta <- paste("provided by user: ", beta) 
    else {
      kb.results$posterior$beta <- list()
      kb.results$posterior$beta$mean <- apply(temp.post$beta.mean, 1, mean)
      names(kb.results$posterior$beta$mean) <- beta.names
      kb.results$posterior$beta$var <- apply(temp.post$beta.var, c(1, 2), mean) + var(t(temp.post$beta.mean))
      dimnames(kb.results$posterior$beta$var) <- list(beta.names,beta.names)
    }
    if(sigmasq.prior == "fixed") kb.results$posterior$sigmasq <- paste("provided by user: ", sigmasq) 
    else{
      kb.results$posterior$sigmasq <- list()
      kb.results$posterior$sigmasq$mean <- mean(temp.post$S2)*df.model/(df.model-2)
      kb.results$posterior$sigmasq$var <- (mean(temp.post$S2)*2/(df.model-4) + var(temp.post$S2))*df.model^2/(df.model-2)^2
    }
    if(phi.prior == "fixed") kb.results$posterior$phi <- paste("provided by user: ", phi) 
    else{
      kb.results$posterior$phi$mean <- mean(kb.results$posterior$phi$sample)
      kb.results$posterior$phi$var <- var(kb.results$posterior$phi$sample)
    }    
    ##
    ## Simulations from the posterior of parameters.
    ##
    if(output$sim.posterior){
      if(beta.size == 1) {
        if(sigmasq.prior == "fixed") {
          if(beta.prior != "fixed")
            kb.results$posterior$beta$sample <- rnorm(n.sim) * as.vector(sqrt(temp.post$beta.var)) + as.vector(temp.post$beta.mean)
        }
        else{
          kb.results$posterior$sigmasq$sample <- rinvchisq(n.sim, df.model, temp.post$S2)
          if(beta.prior != "fixed"){
            cond.beta.sd <- sqrt((as.vector(temp.post$beta.var) * kb.results$posterior$sigmasq$sample)/temp.post$S2)
            kb.results$posterior$beta$sample <- rnorm(n.sim) * cond.beta.sd + as.vector(temp.post$beta.mean)
          }
        }
      }
      else {
        if(sigmasq.prior == "fixed") {
          if(beta.prior != "fixed")
            kb.results$posterior$beta$sample <- array(apply(temp.post$beta.var,3,multgauss),dim=c(beta.size, n.sim))+temp.post$beta.mean
        }
        else {
          kb.results$posterior$sigmasq$sample <- rinvchisq(n.sim, df.model, temp.post$S2)
          if(beta.prior != "fixed"){
            if(is.R()) cond.beta.var <- temp.post$beta.var *rep(kb.results$posterior$sigmasq$sample/temp.post$S2,rep(beta.size^2,n.sim))
            else cond.beta.var <- temp.post$beta.var *rep(kb.results$posterior$sigmasq$sample/temp.post$S2,each = beta.size^2)
            kb.results$posterior$beta$sample <- array(apply(cond.beta.var,3,multgauss),dim=c(beta.size, n.sim)) + temp.post$beta.mean
          }
        }
      }
    }
    if(is.R()) remove("temp.post")
    else remove("temp.post", frame = sys.nframe())
  }
  if(output$keep.mcmc.sim) kb.results$posterior$simulations <- BC.inv(gauss.post, lambda)
  kb.results$model <- model
  kb.results$prior <- prior$priors.info
  kb.results$mcmc.input <- mcmc.input
  kb.results$.Random.seed <- seed
  kb.results$call <- call.fc
  if(is.R()) attr(kb.results, "class") <- c("krige.bayes", "kriging")
  else attr(kb.results, "class") <- setOldClass(c("krige.bayes", "kriging"))
  if(messages.screen) cat("pois.krige.bayes: done!\n")
  return(kb.results)
}
