

"mcmc.binom.aux" <- function(z, data, units.m, meanS, QQ, S.scale, nsim, thin, QtivQ)
{
#
###### ------------------------ doing the mcmc-steps -----------
###### 
#  
  n <- length(data)
  randnormal <- rnorm(n * nsim * thin) * sqrt(S.scale)
  randunif <- runif(nsim * thin)
  z <-  as.double(z)
  S <-  as.double(rep(0, nsim * n))
  acc.rate <-  as.double(1)
  result <- .C("mcmc1binom",
               as.integer(n),
               z = z,
               S = S,
               as.double(data),
	       as.double(units.m),
               as.double(meanS),
               as.double(as.vector(t(QQ))),
               as.double(as.vector(QtivQ)),
               as.double(randnormal),
               as.double(randunif),
               as.double(S.scale),
               as.integer(nsim),
               as.integer(thin),
               acc.rate = acc.rate, DUP=FALSE, PACKAGE = "geoRglm")[c("z", "S", "acc.rate")]
  attr(result$S, "dim") <- c(n, nsim)
  return(result)
}

"mcmc.binom.logit" <- function(data, units.m, meanS, invcov, mcmc.input, messages.screen)
{
####
  ## This is the MCMC engine for the spatial Binomial logit Normal model ----
  ##
  n <- length(data)
  S.scale <- mcmc.input$S.scale
  fisher.l <- data*(1-data/units.m)
  QQ <- t(chol(solve(invcov + diag(fisher.l)))) 
  sqrtfiQ <- sqrt(fisher.l)*QQ 
  QtivQ <- diag(n)-crossprod(sqrtfiQ)  
  if(any(mcmc.input$S.start=="default")) {
    signum <- round(2*data/units.m) -1
    S <- as.vector(ifelse(data > 0 & data < units.m, logit.fct(data/units.m), signum*1.96) - meanS )
    z <- as.vector(solve(QQ,S))
  }
  else{
    if(any(mcmc.input$S.start=="random")) z <- rnorm(n)
    else{
      if(is.numeric(mcmc.input$S.start)){
        if(length(mcmc.input$S.start) != n) stop("dimension of mcmc-starting-value must equal dimension of data")
        else z <- as.vector(solve(QQ,mcmc.input$S.start))
      }
      else stop(" S.start must be a vector of same dimension as data ")
    }
  }
  burn.in <- mcmc.input$burn.in
  thin <- mcmc.input$thin
  n.iter <- mcmc.input$n.iter
## ---------------- burn-in ----------------- ######### 
  if(burn.in > 0) {
    mcmc.output <- mcmc.binom.aux(z, data, units.m, meanS, QQ, S.scale, 1, burn.in, QtivQ)
    if(messages.screen) cat(paste("burn-in = ", burn.in, " is finished. Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate.burn.in <- c(burn.in, mcmc.output$acc.rate)
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
  acc.rate <- matrix(NA, n.turn, 2)
  for(i in 1:n.turn) {
    mcmc.output <- mcmc.binom.aux(mcmc.output$z, data, units.m, meanS, QQ, S.scale, n.temp, thin, QtivQ)
    Sdata[, (n.temp * (i - 1) + 1):(n.temp * i)] <- mcmc.output$S+meanS
    if(messages.screen) cat(paste("iter. numb.", i * n.temp * thin+burn.in, " : Acc.-rate = ", round(mcmc.output$acc.rate, digits=3), "\n"))
    acc.rate[i,1] <-  i * n.temp * thin
    acc.rate[i,2] <- mcmc.output$acc.rate
  }
  if(messages.screen) cat(paste("MCMC performed: n.iter. = ", n.iter, "; thinning = ", thin, "; burn.in = ", burn.in, "\n"))
  if(burn.in > 0) acc.rate <- as.data.frame(rbind(acc.rate.burn.in,acc.rate))
  else acc.rate <- as.data.frame(acc.rate)
  names(acc.rate) <- c("iter.numb", " Acc.rate")
#########
  return(list(Sdata=Sdata, acc.rate=acc.rate))
}


"binom.krige" <- function(geodata, coords = geodata$coords, data = geodata$data, units.m = "default", locations = NULL, borders = NULL, mcmc.input, krige, output)
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
  krige <- krige.glm.check.aux(krige,fct="binom.krige")
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  trend.d <- krige$trend.d
  trend.l <- krige$trend.l
  dist.epsilon <- krige$dist.epsilon
  if(krige$type.krige == "ok") beta.prior <- "flat"
  if(krige$type.krige == "sk") beta.prior <- "deg"
  if(missing(output)) output <- output.glm.control()
  else output <- output.glm.check.aux(output, fct = "binom.krige")
  sim.predict <- output$sim.predict
  messages.screen <- output$messages.screen
  ##
  if(is.vector(coords)){
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
  trend.data <- unclass(trend.spatial(trend=trend.d, geodata = geodata))
  beta.size <- ncol(trend.data)
  if(nrow(trend.data) != n) stop("length of trend is different from the length of the data")
  if(beta.prior == "deg")
    if(beta.size != length(beta))
      stop("size of mean vector is incompatible with trend specified") 
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ## preparing for MCMC 
  ##
  if(missing(mcmc.input)) stop("binom.krige: argument mcmc.input must be given")
  mcmc.input <- mcmc.check.aux(mcmc.input, fct="binom.krige")
  ##
  if(beta.prior == "deg") mean.d <-  as.vector(trend.data %*% beta)
  else mean.d <- rep(0,n)
  if(!is.null(aniso.pars)) {
    invcov <- varcov.spatial(coords = coords.aniso(coords = coords, aniso.pars = aniso.pars), cov.model = cov.model, kappa = kappa, 
                             nugget = nugget, cov.pars = cov.pars, inv = TRUE, func.inv = "cholesky",
                             try.another.decomposition = FALSE)$inverse
  }
  else {
    invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget, cov.pars = cov.pars,
                             inv = TRUE, func.inv = "cholesky", try.another.decomposition = FALSE)$inverse
  }
  ##
########################----- MCMC ------#####################
  ##
  if(beta.prior == "flat") {
    ivtt <- invcov%*%trend.data    
    ittivtt <- solve.geoR(crossprod(trend.data, ivtt))
    invcov <- invcov-ivtt%*%ittivtt%*%t(ivtt)
  }
  prevalence <- mcmc.binom.logit(data = data, units.m = units.m, meanS= mean.d, invcov=invcov, mcmc.input = mcmc.input, messages.screen=messages.screen)
  acc.rate <- prevalence$acc.rate
  prevalence <- plogis(prevalence$Sdata)
  ##
  ##------------------------------------------------------------
######################## ---- prediction ----- #####################
  if(!is.null(locations)) {
    if(!is.null(borders)){
      locations <- locations.inside(locations, borders)
      if(nrow(locations) == 0)
        stop(" binom.krige : there are no prediction locations inside the borders")
      if(messages.screen)
        cat(" binom.krige: results will be returned only for prediction locations inside the borders\n")
    }
    krige <- list(type.krige = krige$type.krige, beta = beta, trend.d = trend.d, trend.l = trend.l, cov.model = cov.model, 
                  cov.pars = cov.pars, kappa = kappa, nugget = nugget, micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
                  aniso.pars = aniso.pars, link = "logit")
    kpl.result <- glm.krige.aux(data = prevalence, coords = coords, locations = locations, krige = krige,
					output = list(n.predictive = ifelse(sim.predict,1,0),
					      signal = TRUE, messages.screen = FALSE))			   
    remove(list = c("prevalence"))
    kpl.result$krige.var <- rowMeans(kpl.result$krige.var) + apply(kpl.result$predict, 1, var) 
    kpl.result$mcmc.error <- sqrt(asympvar(kpl.result$predict)/ncol(kpl.result$predict))
    kpl.result$predict <- rowMeans(kpl.result$predict)
    if(beta.prior == "flat") {
      kpl.result$beta.est <- rowMeans(kpl.result$beta)
      names(kpl.result$beta.est) <- beta.names
    }
    kpl.result$beta <- NULL
  }
  else{
    if(beta.prior == "flat") {
      beta.est <- (ittivtt %*% t(ivtt)) %*% logit.fct(prevalence)
      kpl.result <- list(prevalence=prevalence, beta.est = rowMeans(beta.est), acc.rate=acc.rate)
    }
    else kpl.result <- list(prevalence=prevalence, acc.rate=acc.rate)
  }
  kpl.result$call <- call.fc
#######################################
  attr(kpl.result, "prediction.locations") <- call.fc$locations
  if(!is.null(call.fc$borders)) attr(kpl.result, "borders") <- call.fc$borders
  class(kpl.result) <- "kriging"
  ##class(kpl.result) <- "binom.kriging"
  return(kpl.result)
}

"glm.krige.aux" <- 
function(coords, data, locations, krige, output)
{
  krige$lambda <- 1
  krige$link <- NULL 
  kc.result <- krige.conv.extnd(data = data, coords = coords, locations = locations, krige = krige, output = output)		
  ##
  ##################### Back-transforming predictions
  ##
  predict.transf <- kc.result$predict
  ## using second order taylor-expansion + facts for N(0,1) [third moment = 0 ; fourth moment = 12].
  if(output$messages.screen) cat("binom.krige: back-transforming predictions using 2. order Taylor expansion for g^{-1}() \n")
  ivlogit <- plogis(predict.transf)
  ivlogit1 <- exp(predict.transf)/(1+exp(predict.transf))^2
  ivlogit2 <- exp(predict.transf)*(1-exp(predict.transf))/(1+exp(predict.transf))^3
  ivlogit1[predict.transf>700] <- 0
  ivlogit2[predict.transf>700] <- 0
  kc.result$predict <- ivlogit + 0.5*ivlogit2*kc.result$krige.var
  kc.result$krige.var <- ivlogit1^2*kc.result$krige.var + (11/4)*ivlogit2^2*kc.result$krige.var^2
  remove("predict.transf")
  if(output$n.predictive > 0) {
    kc.result$simulations <- plogis(kc.result$simulations)
  }
  return(kc.result)
}
