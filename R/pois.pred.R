
"mcmc.aux" <- 
  function(z, data, meanS, QQ, Htrunc, S.scale, nsim, thin, Dmat)
{
#
###### ------------------------ doing the mcmc-steps ----------- ############# 
#  
  n <- length(data)
  randnormal <- rnorm(n * nsim * thin) * sqrt(S.scale)
  randunif <- runif(nsim * thin)
  z <-  as.double(z)
  S <-  as.double(rep(0, nsim * n))
  acc.rate <-  as.double(1)
  if(is.null(Dmat)) {
    if(is.R()){
      result <- .C("mcmcrun",
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   acc.rate = acc.rate, DUP=FALSE)[c("z", "S", "acc.rate")]
    }
    else{
      result <- .C("mcmcrun",
                   COPY = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   acc.rate = acc.rate)[c("z", "S", "acc.rate")]
    }
  }
  else {
    if(is.R()){
      result <- .C("mcmcrun2",
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(as.vector(Dmat)),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   acc.rate = acc.rate, DUP=FALSE)[c("z", "S", "acc.rate")]
    }
    else{
      result <- .C("mcmcrun2",
                   COPY = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(as.vector(Dmat)),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   acc.rate = acc.rate)[c("z", "S", "acc.rate")]
    }
  }
  attr(result$S, "dim") <- c(n, nsim)
  return(result)
}


"mcmc.pois.log" <- 
  function(data, units.m, meanS, QQ, mcmc.input, Dmat = NULL)
{
####
  ## This is the MCMC engine for the spatial Poisson log Normal model ----
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n)
      Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  if(any(mcmc.input$S.start=="default")) {
    z <- as.vector(solve(QQ,ifelse(data > 0, log(data), 0) - meanS - log(units.m)))
  }
  else z <- as.vector(solve(QQ,mcmc.input$S.start))
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  ## ---------------- burn-in ----------------- ######### 
  if(burn.in > 0) {
    mcmc.output <- mcmc.aux(z, data, meanS + log(units.m), QQ, Htrunc, S.scale, 1, burn.in, Dmat)
    cat(paste("burn-in = ", burn.in, " is finished. Acc.-rate = ", mcmc.output$acc.rate, "\n"))
  }
  else mcmc.output <- list(z = z)
##### ---------- sampling periode ----------- ###### 
  if(n.iter <= 1000) {
    n.temp <- round(n.iter/thin)
    n.turn <- 1
  }
  else {
    if(thin <= 1000){
      n.temp <- round(1000/thin)
      n.turn <- round(n.iter/1000)
    }
    else{
      n.temp <- 1
      n.turn <- round(n.iter/1000)
    }
  }
  nsim <- n.turn * n.temp
  Sdata <- matrix(NA, n, nsim)
  for(i in 1:n.turn) {
    mcmc.output <- mcmc.aux(mcmc.output$z, data, meanS + log(units.m), QQ, Htrunc, S.scale, n.temp, thin, Dmat)
    Sdata[, (n.temp * (i - 1) + 1):(n.temp * i)] <- mcmc.output$S+meanS
    cat(paste("iter. numb.", i * n.temp * thin, " : Acc.-rate = ", mcmc.output$acc.rate, "\n"))
  }
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  if(is.R()) remove("z")
  else remove(list = c("z"), frame = sys.nframe())
#########
  return(Sdata)
}

"mcmc.boxcox.aux" <- 
  function(z, data, units.m, meanS, QQ, Htrunc, S.scale, nsim, thin, Dmat, lambda)
{
  ##
###### ------------------------ doing the mcmc-steps ----------- ############# 
  ##
  n <- length(data)
  randnormal <- rnorm(n * nsim * thin) * sqrt(S.scale)
  randunif <- runif(nsim * thin)
  z <-  as.double(z)
  S <-  as.double(rep(0, nsim * n))
  acc.rate <-  as.double(1) 
  if(is.null(Dmat)) {
    if(is.R()){
      result <- .C("mcmcrunboxcox",
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(units.m),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   as.double(lambda),
                   acc.rate = acc.rate, DUP=FALSE)[c("z", "S", "acc.rate")]
    }
    else{
      result <- .C("mcmcrunboxcox",
                   COPY = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(units.m),
                   as.double(meanS),
                   as.double(as.vector(t(QQ))),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),  
                   as.double(lambda),
                   acc.rate = acc.rate)[c("z", "S", "acc.rate")]
    }
  }
  else{
    if(is.R()){
      result <- .C("mcmcrun2boxcox",
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(units.m),
                   as.double(as.vector(t(QQ))),
                   as.double(as.vector(Dmat)),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   as.double(lambda),
                   acc.rate = acc.rate, DUP=FALSE)[c("z", "S", "acc.rate")]
    }
    else{
      result <- .C("mcmcrun2boxcox",
                   COPY = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE),
                   as.integer(n),
                   z = z,
                   S = S,
                   as.double(data),
                   as.double(units.m),
                   as.double(as.vector(t(QQ))),
                   as.double(as.vector(Dmat)),
                   as.double(randnormal),
                   as.double(randunif),
                   as.double(Htrunc),
                   as.double(S.scale),
                   as.integer(nsim),
                   as.integer(thin),
                   as.double(lambda),
                   acc.rate = acc.rate)[c("z", "S", "acc.rate")]
    }
  }
  attr(result$S, "dim") <- c(n, nsim)
  return(result)
}

"mcmc.pois.boxcox" <- 
  function(data, units.m, meanS = NULL, QQ, mcmc.input, lambda, Dmat = NULL)
{
####
  ## This is the MCMC engine for the spatial Poisson - Normal model with link from the box-cox-family ----
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  if(is.null(meanS)) meanS <- rep(0,n)
  if(any(mcmc.input$S.start=="default")) {
    S <- as.vector(ifelse(data > 0, ((data/units.m)^lambda-1)/lambda, 0) - meanS )         
    z <- as.vector(solve(QQ,S))
  }
  else z <- as.vector(solve(QQ,mcmc.input$S.start))
  if(any(mcmc.input$Htrunc=="default")) Htrunc <- 2*data + 5
  else {
    if(is.vector(mcmc.input$Htrunc) & length(mcmc.input$Htrunc) == n)
      Htrunc <- mcmc.input$Htrunc
    else Htrunc <- rep(mcmc.input$Htrunc, n)
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
  ## ---------------- burn-in ----------------- ######### 
  if(burn.in > 0) {
    mcmc.output <- mcmc.boxcox.aux(z, data, units.m, meanS, QQ, Htrunc, S.scale, 1, burn.in, Dmat, lambda)
    cat(paste("burn-in = ", burn.in, " is finished. Acc.-rate = ", mcmc.output$acc.rate, "\n"))
  }
  else mcmc.output <- list(z = z)
##### ---------- sampling periode ----------- ###### 
  if(n.iter <= 1000) {
    n.temp <- round(n.iter/thin)
    n.turn <- 1
  }
  else {
    n.temp <- round(1000/thin)
    n.turn <- round(n.iter/1000)
  }
  n.sim <- n.turn * n.temp
  Sdata <- matrix(NA, n, n.sim)
  for(i in 1:n.turn) {
    mcmc.output <- mcmc.boxcox.aux(mcmc.output$z, data, units.m, meanS, QQ, Htrunc, S.scale, n.temp, thin, Dmat, lambda)
    Sdata[, (n.temp * (i - 1) + 1):(n.temp * i)] <- mcmc.output$S+meanS                       
    cat(paste("iter. numb.", i * n.temp * thin, " : Acc.-rate = ", mcmc.output$acc.rate, "\n"))
  }
  cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  if(is.R()) remove("z")
  else remove(list = c("z"), frame = sys.nframe())
#########
  return(Sdata)
}

"krige.glm.control" <-
  function (type.krige = "ok", trend.d = "cte", trend.l = "cte", obj.model = NULL, beta, cov.model, cov.pars, kappa,
            nugget, micro.scale = 0, dist.epsilon = 1e-10, aniso.pars, lambda)
{
  if(type.krige != "ok" & type.krige != "OK" & type.krige != "o.k." & type.krige != "O.K." & type.krige != "sk" & type.krige != "SK" & type.krige != "s.k." & type.krige != "S.K.")
    stop("pois.log.krige: wrong option in the argument type.krige. It should be \"sk\" or \"ok\"(if ordinary or simple kriging is to be performed)")
  if(type.krige=="OK" | type.krige=="O.K." |type.krige=="o.k.")
    type.krige <- "ok"
  if(type.krige=="SK" | type.krige=="S.K." |type.krige=="s.k.")
    type.krige <- "sk"
  ##
  if(!is.null(obj.model)){
    if(missing(beta)) beta <- obj.model$beta
    if(missing(cov.model)) cov.model <- obj.model$cov.model
    if(missing(cov.pars)) cov.pars <- obj.model$cov.pars
    if(missing(kappa)) kappa <- obj.model$kappa
    if(missing(nugget)) nugget <- obj.model$nugget
    if(missing(lambda)) lambda <- obj.model$lambda
    if(missing(aniso.pars)) aniso.pars <- obj.model$aniso.pars
  }
  else{
    if(missing(beta)) beta <- NULL
    if(missing(cov.model)) cov.model <- "matern"
    if(missing(cov.pars))
      stop("covariance parameters (sigmasq and phi) should be provided")
    if(missing(kappa)) kappa <- 0.5
    if(missing(nugget)) nugget <- 0
    if(missing(lambda)) lambda <- 0
    if(missing(aniso.pars)) aniso.pars <- NULL
  }
  ##
  if(type.krige == "sk")
    if(is.null(beta) | !is.numeric(beta))
      stop("\n pois.log.krige: argument beta must be provided in order to perform simple kriging")
  if(micro.scale > nugget)
    stop("pois.log.krige: micro.scale must be in the interval [0, nugget]")
  if(!is.null(aniso.pars))
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("pois.log.krige: anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
  ##
  if(inherits(trend.d, "formula") | inherits(trend.l, "formula")){
    if((!inherits(trend.d, "formula")) | (!inherits(trend.l, "formula")))
      stop("pois.log.krige: trend.d and trend.l must have similar specification")
  }
  else{
    if(trend.d != trend.l)
      stop("pois.log.krige: trend.l is different from trend.d")
  }
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential","gaussian",
                           "spherical", "circular", "cubic",
                           "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  res <- list(type.krige = type.krige,
              trend.d = trend.d, trend.l = trend.l, 
              beta = beta,
              cov.model = cov.model, 
              cov.pars = cov.pars, kappa = kappa,
              nugget = nugget,
              micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
              aniso.pars = aniso.pars, lambda = lambda)
  class(res) <- "krige.geoRglm"
  return(res)
}

"pois.log.krige" <- 
function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", locations = NULL,
         mcmc.input, krige, output)
{
  if(missing(geodata))
    geodata <- list(coords=coords, data=data)
  call.fc <- match.call()
  n <- length(data)
  if(any(units.m == "default")){
    if(!is.null(geodata$units.m)) units.m <- geodata$units.m
    else units.m <- rep(1, n)
  }
  if(missing(krige)) stop("must provide object krige")
  else{
    if(is.null(class(krige)) || class(krige) != "krige.geoRglm"){
      if(!is.list(krige))
        stop("pois.log.krige: the argument krige only takes a list or an output of the function krige.glm.control")
      else{
        krige.names <-c("type.krige","trend.d","trend.l","obj.model","beta","cov.model",
"cov.pars","kappa","nugget","micro.scale","dist.epsilon","lambda","aniso.pars")
        krige.user <- krige
        krige <- list()
        if(length(krige.user) > 0){
          for(i in 1:length(krige.user)){
            n.match <- match.arg(names(krige.user)[i], krige.names)
            krige[[n.match]] <- krige.user[[i]]
          }
        }
        if(is.null(krige$type.krige)) krige$type.krige <- "ok"  
        if(is.null(krige$trend.d)) krige$trend.d <-  "cte"
        if(is.null(krige$trend.l)) krige$trend.l <-  "cte"
        if(is.null(krige$cov.model)) krige$cov.model <- "matern"
        if(is.null(krige$kappa)) krige$kappa <-  0.5
        if(is.null(krige$nugget)) krige$nugget <-  0
        if(is.null(krige$micro.scale)) krige$micro.scale <- 0  
        if(is.null(krige$dist.epsilon)) krige$dist.epsilon <-  1e-10
        if(is.null(krige$lambda)) krige$lambda <- 0
          krige <- krige.glm.control(type.krige = krige$type.krige,
                                 trend.d = krige$trend.d, trend.l = krige$trend.l,
                                 obj.model = krige$obj.model,
                                 beta = krige$beta, cov.model = krige$cov.model,
                                 cov.pars = krige$cov.pars, kappa = krige$kappa,
                                 nugget = krige$nugget, micro.scale = krige$micro.scale,
                                 dist.epsilon = krige$dist.epsilon, 
                                 aniso.pars = krige$aniso.pars,
                                 lambda = krige$lambda)
        
      }
    }
  }
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  if(krige$lambda>0.001) stop("pois.log.krige : lambda > 0 is not implemented ")
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  if(krige$micro.scale>0.00000001) stop("pois.log.krige : microscale > 0 is not implemented ")
  aniso.pars <- krige$aniso.pars
  trend.d <- krige$trend.d
  trend.l <- krige$trend.l
  dist.epsilon <- krige$dist.epsilon
  if(krige$type.krige == "ok") beta.prior <- "flat"
  if(krige$type.krige == "sk") beta.prior <- "deg"
  if(missing(output))
    output <- output.glm.control()
  else{
    if(is.null(class(output)) || class(output) != "output.geoRglm"){
      if(!is.list(output))
        stop("pois.log.krige: the argument output only takes a list or an output of the function output.glm.control")
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
        if(is.null(output$sim.predict)) output$sim.predict <- FALSE      
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE
        output <- output.glm.control(sim.predict = output$sim.predict,
                                 messages.screen = output$messages.screen)
      }
    }
  }
  sim.predict <- output$sim.predict
  messages.screen <- output$messages.screen
  ##
  if(is.vector(coords)) {
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimension assumed")
  }
  coords <- as.matrix(coords)
  dimnames(coords) <- list(NULL, NULL)
  ##
  if(is.null(locations)) {
    cat(paste("locations need to be specified for prediction; prediction not performed \n"))
  }
  else {
    if(is.null(trend.l))
      stop("trend.l needed for prediction")
  }
  trend.data <- trend.spatial(trend = trend.d, geodata=geodata)
  beta.size <- ncol(trend.data)
  if(beta.prior == "deg")
    if(beta.size != length(beta))
      stop("size of mean vector is incompatible with trend specified") 
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ## preparing for MCMC 
  ##
  if(missing(mcmc.input)) stop("pois.log.krige: argument mcmc.input must be given")
  else{
    if(is.null(class(mcmc.input)) || class(mcmc.input) != "mcmc.geoRglm"){
      if(!is.list(mcmc.input))
        stop("pois.log.krige: the argument mcmc.input only takes a list or an output of the function mcmc.control")
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
        mcmc <- mcmc.control(S.scale = mcmc.input$S.scale,Htrunc=mcmc.input$Htrunc,S.start=mcmc.input$S.start,
                             burn.in=mcmc.input$burn.in,thin=mcmc.input$thin,n.iter=mcmc.input$n.iter,
                             phi.start=mcmc.input$phi.start,phi.scale=mcmc.input$phi.scale)
      }
    }
  }
  ##
  if(beta.prior == "deg") mean.d <-  as.vector(trend.data %*% beta)
  else mean.d <- 0
  if(!is.null(aniso.pars)) {
    QQt <- varcov.spatial(coords = coords.aniso(coords = coords, aniso.pars = aniso.pars), cov.model = cov.model, kappa = kappa, 
           nugget = nugget, cov.pars = cov.pars, only.decomposition = TRUE, func.inv = "cholesky", 
           try.another.decomposition = FALSE)$sqrt.varcov
  }
  else {
    QQt <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget, cov.pars = cov.pars,
            only.decomposition = TRUE, func.inv = "cholesky", try.another.decomposition = FALSE)$sqrt.varcov
  }
  ##
########################----- MCMC ------#####################
  ##
  if(beta.prior == "deg") {
    intensity <- exp(mcmc.pois.log(data = data, units.m = units.m, meanS = mean.d, QQ = t(QQt), mcmc.input = mcmc.input))
    if(is.R()) remove(list = c("QQt"))
    else remove(list = c("QQt"), frame = sys.nframe())
  }
  if(beta.prior == "flat") {
    if(is.R()){
      QQtinv <- solve(QQt)
      invcov <- chol2inv(QQt)
    }
    else{
      QQtinv <- solve.upper(QQt)
      invcov <- QQtinv %*% t(QQtinv)
    }
    temp <- crossprod(QQtinv,trend.data)
    ittivtt <- solve.geoR(crossprod(temp))
    intensity <- exp(mcmc.pois.log(data = data, units.m = units.m, meanS = mean.d,QQ = t(QQt), mcmc.input = mcmc.input,
                                   Dmat = diag(n) - temp %*% ittivtt %*% t(temp)))
    if(is.R()) remove(list = c("QQt", "QQtinv", "temp"))
    else
      remove(list = c("QQt", "QQtinv", "temp"), frame = sys.nframe())
  }
  if(is.R()) remove(list = c("mean.d"))
  else remove(list = c("mean.d"), frame = sys.nframe())
  ##
  ##------------------------------------------------------------
######################## ---- prediction ----- #####################
  if(!is.null(locations)) {
    krige <- list(type.krige = krige$type.krige, beta = beta, trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, 
             cov.pars = cov.pars, kappa = kappa, nugget = nugget, micro.scale = nugget, dist.epsilon = dist.epsilon, 
             aniso.pars = aniso.pars, lambda = 0)
    kpl.result <- krige.conv.extnd(data = intensity, coords = coords, locations = locations, krige = krige,
                                   output = list(n.predictive = ifelse(sim.predict,1,0), signal = FALSE, messages.screen = FALSE))
    if(is.R()) remove(list = c("intensity"))
    else remove(list = c("intensity"), frame = sys.nframe())
    kpl.result$krige.var <- apply(kpl.result$krige.var, 1, mean) + apply(kpl.result$predict, 1, var)
    kpl.result$mcmc.error <- sqrt(asympvar(kpl.result$predict)/ncol(kpl.result$predict))
    kpl.result$predict <- apply(kpl.result$predict, 1, mean)
    ##if(!krige$signal)
      ##kpl.result$krige.var <- kpl.result$krige.var + kpl.result$predict    #### to be done properly sometime !
    if(beta.prior == "flat") {
      kpl.result$beta.est <- apply(kpl.result$beta, 1, mean)
      names(kpl.result$beta.est) <- beta.names
    }
    kpl.result$beta <- NULL
  }
  else{
    if(beta.prior == "flat") {
      beta.est <- (ittivtt %*% crossprod(trend.data, invcov)) %*% log(intensity)
      kpl.result <- list(intensity=intensity, beta.est = apply(beta.est, 1, mean))
    }
    else kpl.result <- list(intensity=intensity)
  }
  kpl.result$call <- call.fc
#######################################
  attr(kpl.result, "prediction.locations") <- call.fc$locations
  class(kpl.result) <- "kriging"
  return(kpl.result)
}
