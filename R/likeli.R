"prepare.lik.sim" <-
  function(bayes.output,geodata, use.intensity = FALSE)
{
  n.dat <- nrow(bayes.output$posterior$simulations)
  n.sim <- ncol(bayes.output$posterior$simulations)
  if(use.intensity){
    if(any(geodata$data == 0)) stop("use.intensity = TRUE is only allowed when all data are positive ")
  }
  lambda <- bayes.output$model$lambda
  if(lambda == 0)
    S <- log(bayes.output$posterior$simulations)
  else 
    S <- (bayes.output$posterior$simulations^lambda-1)/lambda
  if(bayes.output$model$trend.d == "cte") trend <- rep(1,n.dat)
  else{
    if(is.matrix(bayes.output$model$trend.d)) stop("only supporting one-dimensional beta")
    else trend <- bayes.output$model$trend.d
  }
  beta.size <- 1
  cov.model <- bayes.output$model$cov.model
  kappa <- bayes.output$model$kappa
  nugget.rel <- bayes.output$prior$tausq.rel$fixed.value
  if(is.null(bayes.output$prior$beta$dist) || bayes.output$prior$beta$dist != "flat") stop("Only object with flat beta prior is allowed")
  if(!is.null(bayes.output$prior$sigmasq$dist) && bayes.output$prior$sigmasq$dist == "uniform"){
    n.sigma <- -2
    s2.sigma <- 0
  }
  else stop("OFC, JEG ER IKKE SIKKER PAA AT DET HER VIRKER ")
  if(is.null(bayes.output$prior$phi$status) || bayes.output$prior$phi$status != "fixed") stop("Only object with fixed phi is allowed")
  else phi <- bayes.output$prior$phi$fixed.value
  if(is.null(bayes.output$model$aniso.pars)) coords <- geodata$coords
  else coords <- coords.aniso(coords = geodata$coords, aniso.pars = bayes.output$model$aniso.pars)
  log.f.sim <- rep(0,n.sim)
  invcov <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget.rel, 
                           cov.pars = c(1,phi), det = TRUE, only.inv.lower.diag = TRUE)
  temp3 <- diagquadraticformXAX(t(trend), invcov$lower.inverse, invcov$diag.inverse)  
  temp1 <- diagquadraticformXAX(S, invcov$lower.inverse, invcov$diag.inverse)
  temp2 <- as.vector(bilinearformXAY(S, invcov$lower.inverse, invcov$diag.inverse, trend)) 
  log.f.sim <- 0.5*log(temp3)-invcov$log.det.to.half - 0.5*(n.dat-beta.size+n.sigma)*log(n.sigma*s2.sigma+temp1-temp2^2/temp3)
  if(use.intensity){
    logJ <- apply(log(bayes.output$posterior$simulations),2,sum)*(lambda-1)
    return(list(mu = bayes.output$posterior$simulations, coords = geodata$coords, aniso.pars = bayes.output$model$aniso.pars,
                lambda = lambda, log.f.sim = log.f.sim + logJ))
  }
  else{
    return(list(S = S, coords = coords, aniso.pars = bayes.output$model$aniso.pars, lambda = lambda, log.f.sim = log.f.sim))
  }  
}

"maxim.aux1" <-
function(S,invcov,trend = "cte",log.f.sim)
{
   n.dat <- nrow(S) 
   if(trend == "cte") trend <- rep(1, n.dat)
   if(is.matrix(trend)) stop("only supporting one-dimensional beta") 
   temp1 <- diagquadraticformXAX(S,invcov$lower.inverse,invcov$diag.inverse)
   temp2 <- as.vector(bilinearformXAY(S,invcov$lower.inverse,invcov$diag.inverse,trend))
   temp3 <- as.vector(bilinearformXAY(trend,invcov$lower.inverse,invcov$diag.inverse,trend))
   beta.hat <- temp2/temp3
   sigmasq.hat <- diagquadraticformXAX(S-trend%o%beta.hat,invcov$lower.inverse,invcov$diag.inverse)/n.dat
   beta <- beta.hat[1]
   sigmasq <- sigmasq.hat[1]
   corr1 <- mean(0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq-log.f.sim)
   log.f.sim <- log.f.sim + corr1
   hh <- 0 
   for(ll in 1:ncol(S)){
      functemp <- exp(-0.5*(temp1-2*temp2*beta.hat[ll]+temp3*beta.hat[ll]^2)/sigmasq.hat[ll]-log.f.sim)
      hhnew <- mean(functemp)/(sigmasq.hat[ll]^(n.dat/2))
      if(hhnew>hh){
         beta <-beta.hat[ll]
         sigmasq <-sigmasq.hat[ll]
         hh <- hhnew   
      }
   }
   functemp <- exp(-0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq-log.f.sim)
   test <- 1
   test2 <- 1
   steplen <- 1
   while(test>0.0000000000001 | test2 > 0 ){
      F1 <- mean((temp2-temp3*beta)*functemp)/(sigmasq^(n.dat/2+1)) 
      F2 <- -(n.dat/2)*hh/sigmasq + 0.5*mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
      F11 <- mean((temp2-temp3*beta)^2*functemp)/(sigmasq^(n.dat/2+2))-hh*temp3/sigmasq 
      F12 <- -(n.dat/2+1)*F1/sigmasq+0.5*mean((temp2-temp3*beta)*(temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
      F22 <- (n.dat/2)*(n.dat/2+1)*hh/(sigmasq^(n.dat/2+2))-mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+3)) + 0.25*mean((temp1-2*temp2*beta+temp3*beta^2)^2*functemp)/(sigmasq^(n.dat/2+4))
      Deter <-  (F11*F22-F12*F12)      
      if(Deter != 0){
        betanew <- beta - steplen*(F22*F1 - F12*F2)/Deter
        sigmasqnew <- sigmasq - steplen*(-F12*F1 + F11*F2)/Deter
      }
      else {
        cat(paste("Problems when optimising w.r.t. beta and sigmasq: likelihood is very flat \n"))
        cat(paste("likelihood value at this stage is = ",log(hh)+corr1,"\n"))
        betanew <- beta 
        sigmasqnew <- sigmasq 
      }
      functempnew <- exp(-0.5*(temp1-2*temp2*betanew+temp3*betanew^2)/sigmasqnew-log.f.sim)
      hhnew <- mean(functempnew)/(sigmasqnew^(n.dat/2)) 
      test <- (beta-betanew)^2+(sigmasq-sigmasqnew)^2 
      test2 <- hhnew-hh 
      if(test2>0){ 
         hh <- hhnew  
         beta <- betanew  
         sigmasq <- sigmasqnew  
         functemp <- functempnew  
      }
      else{
         steplen <- steplen/2     
      }
   }
   return(list(beta = beta, sigmasq = sigmasq, logh = log(hh)+corr1)) 
}

"lik.sim" <-
function(pars, fp, ip, temp.list)
{ 
  if(is.R()) require(mva)
  n.dat <- temp.list$n
  beta.size <- temp.list$beta.size
  if(beta.size>1) stop("only supporting one-dimensional beta")
  temp.list$xmat <- as.vector(temp.list$xmat)
  ## Obligatory parameter:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  ##
  S <- temp.list$z
  ##
  ## Computing likelihood
  ##
  iv <- varcov.spatial(dists.lowertri = as.vector(dist(temp.list$coords)), cov.model = temp.list$cov.model, kappa = temp.list$kappa,
                         nugget = tausq.rel, cov.pars = c(1, phi), only.inv.lower.diag = TRUE, det = TRUE)
  temp1 <- diagquadraticformXAX(S,iv$lower.inverse,iv$diag.inverse)
  temp2 <- as.vector(bilinearformXAY(S, iv$lower.inverse, iv$diag.inverse,temp.list$xmat))
  temp3 <- as.vector(bilinearformXAY(temp.list$xmat,iv$lower.inverse,iv$diag.inverse,temp.list$xmat))
  beta.hat <- temp2/temp3 
  sigmasq.hat <- diagquadraticformXAX(S-temp.list$xmat%o%beta.hat,iv$lower.inverse,iv$diag.inverse)/n.dat
  beta <- beta.hat[1]
  sigmasq <- sigmasq.hat[1]
  corr1 <- mean(0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq-temp.list$log.f.sim)
  log.f.sim <- temp.list$log.f.sim + corr1
  hh <- 0 
  if(is.null(temp.list$messages.screen)){
     cat(paste(phi,tausq.rel,"\n"))
     messages.screen <- TRUE
  }
  else{
     messages.screen <- temp.list$messages.screen
     if(messages.screen) cat(paste(phi,tausq.rel,"\n"))
  }  
  for(ll in 1:ncol(S)){
      functemp <- exp(-0.5*(temp1-2*temp2*beta.hat[ll]+temp3*beta.hat[ll]^2)/sigmasq.hat[ll]-log.f.sim)
      hhnew <- mean(functemp)/(sigmasq.hat[ll]^(n.dat/2))
      if(hhnew>hh){
         beta <-beta.hat[ll]
         sigmasq <-sigmasq.hat[ll]
         hh <- hhnew 
      }
  } 
  functemp <- exp(-0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq-log.f.sim)
  test <- 1
  test2 <- 1
  steplen <- 1
  while(test>0.0000000000001 | test2 > 0 ){
      F1 <- mean((temp2-temp3*beta)*functemp)/(sigmasq^(n.dat/2+1)) 
      F2 <- -(n.dat/2)*hh/sigmasq + 0.5*mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
      F11 <- mean((temp2-temp3*beta)^2*functemp)/(sigmasq^(n.dat/2+2))-hh*temp3/sigmasq 
      F12 <- -(n.dat/2+1)*F1/sigmasq+0.5*mean((temp2-temp3*beta)*(temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
      F22 <- (n.dat/2)*(n.dat/2+1)*hh/(sigmasq^(n.dat/2+2))-mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2++3)) + 0.25*mean((temp1-2*temp2*beta+temp3*beta^2)^2*functemp)/(sigmasq^(n.dat/2+4))
      Deter <-  (F11*F22-F12*F12)
      if(Deter != 0){
        betanew <- beta - steplen*(F22*F1 - F12*F2)/Deter
        sigmasqnew <- sigmasq - steplen*(-F12*F1 + F11*F2)/Deter
      }
      else {
        if(messages.screen){
           cat(paste("Problems when optimising w.r.t. beta and sigmasq: likelihood is very flat  \n"))
           cat(paste("likelihood value at this stage is = ",log(hh)+corr1 - iv$log.det.to.half,"\n"))
        }
        betanew <- beta 
        sigmasqnew <- sigmasq 
      }
      functempnew <- exp(-0.5*(temp1-2*temp2*betanew+temp3*betanew^2)/sigmasqnew-log.f.sim)
      hhnew <- mean(functempnew)/(sigmasqnew^(n.dat/2)) 
      test <- (beta-betanew)^2+(sigmasq-sigmasqnew)^2 
      test2 <- hhnew-hh 
      if(test2>0){ 
         hh <- hhnew  
         beta <- betanew  
         sigmasq <- sigmasqnew  
         functemp <- functempnew  
      }
      else{
         steplen <- steplen/2     
      }
  }
  # calculating M.L.E.
  #
  negloglik <- ( - log(hh) - corr1 + iv$log.det.to.half)
  if(messages.screen) cat(paste(-negloglik,"\n"))
  return(negloglik)
}

"lik.sim.boxcox" <-
function(pars, fp, ip, temp.list)
{ 
### Function for finding m.l.e. for a given phi based on samples from mu=g^{-1}(S) ###############
### This function is only valid when all observations are positive.
##
  if(is.R()) require(mva)
  n.dat <- temp.list$n
  beta.size <- temp.list$beta.size
  if(beta.size>1) stop("only supporting one-dimensional beta")
  temp.list$xmat <- as.vector(temp.list$xmat)
  ## Obligatory parameter:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  ##
  mu <- temp.list$mu
  ##
  ## computing the determinant of the transformation
  ##
  if(ip$f.lambda) lambda <- fp$lambda
  else lambda <- pars[length(pars)]
  log.J.lambda <- apply(log(mu),2,sum)*(lambda-1)
  if(lambda ==0) mu <- log(mu)
  else mu <- (mu^lambda-1)/lambda
  ##
  ## Computing likelihood
  ##
  iv <- varcov.spatial(dists.lowertri = as.vector(dist(temp.list$coords)), cov.model = temp.list$cov.model, kappa = temp.list$kappa,
                         nugget = tausq.rel, cov.pars = c(1, phi), only.inv.lower.diag = TRUE, det = TRUE)
  temp1 <- diagquadraticformXAX(mu,iv$lower.inverse,iv$diag.inverse)
  temp2 <- as.vector(bilinearformXAY(mu, iv$lower.inverse, iv$diag.inverse,temp.list$xmat))
  temp3 <- as.vector(bilinearformXAY(temp.list$xmat,iv$lower.inverse,iv$diag.inverse,temp.list$xmat))
  beta.hat <- temp2/temp3 
  sigmasq.hat <- diagquadraticformXAX(mu-temp.list$xmat%o%beta.hat,iv$lower.inverse,iv$diag.inverse)/n.dat
  beta <- beta.hat[1]
  sigmasq <- sigmasq.hat[1]
  corr1 <- mean(0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq+log.J.lambda-temp.list$log.f.sim)
  log.f <- temp.list$log.f.sim - log.J.lambda + corr1
  hh <- 0 
  if(is.null(temp.list$messages.screen)){
     cat(paste(phi,tausq.rel,lambda,"\n"))
     messages.screen <- TRUE
  }
  else{
     if(temp.list$messages.screen) cat(paste(phi,tausq.rel,lambda,"\n"))
     messages.screen <- temp.list$messages.screen
  }
  for(ll in 1:ncol(mu)){
     functemp <- exp(-0.5*(temp1-2*temp2*beta.hat[ll]+temp3*beta.hat[ll]^2)/sigmasq.hat[ll]-log.f)
     hhnew <- mean(functemp)/(sigmasq.hat[ll]^(n.dat/2))
     if(hhnew>hh){
         beta <- beta.hat[ll]
         sigmasq <- sigmasq.hat[ll]
         hh <- hhnew 
     }
  }
  functemp <- exp(-0.5*(temp1-2*temp2*beta+temp3*beta^2)/sigmasq-log.f)
  test <- 1
  test2 <- 1
  steplen <- 1
  while(test>0.0000000000001 | test2 > 0 ){
     F1 <- mean((temp2-temp3*beta)*functemp)/(sigmasq^(n.dat/2+1)) 
     F2 <- -(n.dat/2)*hh/sigmasq + 0.5*mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
     F11 <- mean((temp2-temp3*beta)^2*functemp)/(sigmasq^(n.dat/2+2))-hh*temp3/sigmasq 
     F12 <- -(n.dat/2+1)*F1/sigmasq+0.5*mean((temp2-temp3*beta)*(temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2+2)) 
     F22 <- (n.dat/2)*(n.dat/2+1)*hh/(sigmasq^(n.dat/2+2))-mean((temp1-2*temp2*beta+temp3*beta^2)*functemp)/(sigmasq^(n.dat/2++3)) + 0.25*mean((temp1-2*temp2*beta+temp3*beta^2)^2*functemp)/(sigmasq^(n.dat/2+4))
     Deter <-  (F11*F22-F12*F12)
     if(Deter != 0){
        betanew <- beta - steplen*(F22*F1 - F12*F2)/Deter
        sigmasqnew <- sigmasq - steplen*(-F12*F1 + F11*F2)/Deter
     }
     else {
        if(messages.screen){
           cat(paste("Problems when optimising w.r.t. beta and sigmasq: likelihood is very flat  \n"))
           cat(paste("likelihood value at this stage is = ",hh,"\n"))
        }
        betanew <- beta 
        sigmasqnew <- sigmasq 
     }
     functempnew <- exp(-0.5*(temp1-2*temp2*betanew+temp3*betanew^2)/sigmasqnew-log.f)
     hhnew <- mean(functempnew)/(sigmasqnew^(n.dat/2)) 
     test <- (beta-betanew)^2+(sigmasq-sigmasqnew)^2 
     test2 <- hhnew-hh 
     if(test2>0){ 
         hh <- hhnew  
         beta <- betanew  
         sigmasq <- sigmasqnew  
         functemp <- functempnew  
     }
     else{
        steplen <- steplen/2     
     }
  }
  # calculating M.L.E.
  #
  negloglik <- ( - log(hh) - corr1 + iv$log.det.to.half)
  if(messages.screen) cat(paste(-negloglik,"\n")) 
  return(negloglik)
}

"intgr.lik.sim" <-
function(pars, fp, ip, temp.list)
{ 
  if(is.R()) require(mva)
  n.dat <- temp.list$n
  beta.size <- temp.list$beta.size
  if(beta.size>1) stop("only supporting one-dimensional beta")
  temp.list$xmat <- as.vector(temp.list$xmat)
  ## Obligatory parameter:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  ##
  data <- temp.list$z
  if(is.null(temp.list$messages.screen)){
     cat(paste(phi,tausq.rel,"\n"))
     messages.screen <- TRUE
  }
  else{
     if(temp.list$messages.screen) cat(paste(phi,tausq.rel,"\n"))
     messages.screen <- temp.list$messages.screen
  }
  ##
  ## Computing integrated likelihood
  ##
  n.sigma <- -2
  s2.sigma <- 0
  invcov <- varcov.spatial(coords = temp.list$coords, cov.model = temp.list$cov.model, kappa = temp.list$kappa, nugget = tausq.rel, 
                    cov.pars = c(1,phi), det = TRUE, only.inv.lower.diag = TRUE)
  temp3 <- as.vector(bilinearformXAY(temp.list$xmat,invcov$lower.inverse,invcov$diag.inverse,temp.list$xmat))
  temp1 <- diagquadraticformXAX(data,invcov$lower.inverse,invcov$diag.inverse)
  temp2 <- as.vector(bilinearformXAY(data,invcov$lower.inverse,invcov$diag.inverse,temp.list$xmat)) 
  temp.func <- (n.sigma*s2.sigma+temp1-temp2^2/temp3)^(-0.5*(n.dat-beta.size+n.sigma))*exp(-temp.list$log.f.sim)
  negloglik <- ( - log(mean(temp.func)) + invcov$log.det.to.half - 0.5*log(temp3) )
  if(messages.screen) cat(paste(-negloglik,"\n"))
  return(negloglik)
}

"intgr.lik.sim.boxcox" <-
function(pars, fp, ip, temp.list)
{ 
  if(is.R()) require(mva)
  n.dat <- temp.list$n
  beta.size <- temp.list$beta.size
  if(beta.size>1) stop("only supporting one-dimensional beta")
  temp.list$xmat <- as.vector(temp.list$xmat)
  ## Obligatory parameter:
  phi <- pars[1]
  if(ip$f.tausq.rel) tausq.rel <- fp$tausq.rel
  else tausq.rel <- pars[2]
  ##
  mu <- temp.list$mu
  ##
  ## computing the determinant of the transformation
  ##
  if(ip$f.lambda) lambda <- fp$lambda
  else lambda <- pars[length(pars)]
  log.J.lambda <- apply(log(mu),2,sum)*(lambda-1)
  if(lambda ==0) mu <- log(mu)
  else mu <- (mu^lambda-1)/lambda
  if(is.null(temp.list$messages.screen)){
     cat(paste(phi,tausq.rel,lambda,"\n"))
     messages.screen <- TRUE
  }
  else{
     if(temp.list$messages.screen) cat(paste(phi,tausq.rel,lambda,"\n"))
     messages.screen <- temp.list$messages.screen
  }
  ##
  ## Computing integrated likelihood
  ##
  n.sigma <- -2
  s2.sigma <- 0
  invcov <- varcov.spatial(coords = temp.list$coords, cov.model = temp.list$cov.model, kappa = temp.list$kappa, nugget = tausq.rel, 
                    cov.pars = c(1,phi), det = TRUE, only.inv.lower.diag = TRUE)
  temp3 <- as.vector(bilinearformXAY(temp.list$xmat,invcov$lower.inverse,invcov$diag.inverse,temp.list$xmat))
  temp1 <- diagquadraticformXAX(mu,invcov$lower.inverse,invcov$diag.inverse)
  temp2 <- as.vector(bilinearformXAY(mu,invcov$lower.inverse,invcov$diag.inverse,temp.list$xmat)) 
  temp.func <- (n.sigma*s2.sigma+temp1-temp2^2/temp3)^(-0.5*(n.dat-beta.size+n.sigma))*exp(-temp.list$log.f.sim + log.J.lambda)
  negloglik <- ( - log(mean(temp.func)) + invcov$log.det.to.half - 0.5*log(temp3))
  if(messages.screen) cat(paste(-negloglik,"\n"))
  return(negloglik)
}

"pois.likfit.sim" <-
function (mcmc.obj, trend = "cte", 
    cov.model = c("matern", "exponential", "gaussian","spherical", "circular", "cubic", "wave",
            "powered.exponential", "cauchy", "gneiting","gneiting.matern", "pure.nugget"), 
    kappa = 0.5, ini.phi, fix.nugget.rel = FALSE, nugget.rel = 0, aniso.pars = NULL, 
    fix.lambda = TRUE, lambda = NULL, lik.type = "standard", limits = likfit.limits(), messages.screen = TRUE, ...)
{
  ##
  ## Checking input
  ##
    geodata <- list(coords=mcmc.obj$coords)
  if(is.R()) require(mva)
  call.fc <- match.call()
  temp.list <- list()
  ##
  cov.model <- match.arg(cov.model)
  if(!is.null(kappa)){
    if(cov.model == "matern" & kappa == 0.5) cov.model <- "exponential"
  }
  ##
  if(is.null(mcmc.obj$S)){
     if(is.null(mcmc.obj$mu)) stop("mcmc.obj should include either an object mu or an object S.") 
     if(is.null(lambda)) stop("must specify lambda")
     n <- temp.list$n <- nrow(mcmc.obj$mu)
     temp.list$mu <- mcmc.obj$mu
     boxcox <- TRUE
  }
  else{
     if(!is.null(lambda)) cat(paste("cannot use argument lambda with the given objects in mcmc.obj"))
     if(!fix.lambda){
        cat(paste("cannot estimate lambda from the given objects in mcmc.obj"))
        fix.lambda <- TRUE
     }
     n <- temp.list$n <- nrow(mcmc.obj$S)
     temp.list$z <- mcmc.obj$S 
     boxcox <- FALSE
  }
  if (n != nrow(mcmc.obj$coords)) stop("Number of locations does not match with number of data")
  if(!is.null(aniso.pars)){
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
       stop("anisotropy parameters must be a vector with two elements: rotation angle (in radians) and anisotropy ratio (a number > 1)")
    coords <- coords.aniso(coords = as.matrix(mcmc.obj$coords), aniso.pars = aniso.pars)
  }
  else{
  
     if(!is.null(mcmc.obj$aniso.pars)){
        coords <- coords.aniso(coords = as.matrix(mcmc.obj$coords), aniso.pars = mcmc.obj$aniso.pars)
        aniso.pars <- mcmc.obj$aniso.pars
     }
     else coords <- as.matrix(mcmc.obj$coords)
  }
  temp.list$xmat <- trend.spatial(trend = trend, geodata=geodata)
  beta.size <- temp.list$beta.size <- dim(temp.list$xmat)[2]  
  ##
  temp.list$coords <- coords
  temp.list$cov.model <- cov.model
  temp.list$kappa <- kappa
  ini <- ini.phi
  lower.optim <- c(limits$phi["lower"])
  upper.optim <- c(limits$phi["upper"])
  fixed.values <- list()
  ##
  if(fix.nugget.rel) {
    fixed.values$tausq.rel <- nugget.rel
  }
  else {
    ini <- c(ini, nugget.rel)
    lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
    upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
  }
  if(boxcox){
     if(fix.lambda) {
       fixed.values$lambda <- lambda
     }
     else {
       ini <- c(ini, lambda)
       lower.optim <- c(lower.optim, limits$lambda["lower"])
       upper.optim <- c(upper.optim, limits$lambda["upper"])
     }
     ip <- list(f.tausq.rel = fix.nugget.rel, f.lambda = fix.lambda)
  }
  else ip <- list(f.tausq.rel = fix.nugget.rel)
  names(ini) <- NULL
  temp.list$log.f.sim <- mcmc.obj$log.f.sim
  temp.list$messages.screen <- messages.screen
  ##
  if(messages.screen){
    cat("--------------------------------------------------------------------\n")
    cat("pois.likfit.sim: likelihood maximisation using the function ")
    if(is.R()) cat("optim.\n") else cat("nlminb.\n")     
  }
  npars <- beta.size + 2 + sum(!unlist(ip))
  if(lik.type == "standard"){
     if(is.R()){
        if(boxcox){
          lik.optim <- optim(par = ini, fn = lik.sim.boxcox, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
        else{
           lik.optim <- optim(par = ini, fn = lik.sim, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
     }
     else{
        if(boxcox){
           lik.optim <- nlminb(ini, lik.sim.boxcox, lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
        else{
           lik.optim <- nlminb(ini, lik.sim, lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
     }
  }
  else{
     if(is.R()){
        if(boxcox){
          lik.optim <- optim(par = ini, fn = intgr.lik.sim.boxcox, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
        else{
           lik.optim <- optim(par = ini, fn = intgr.lik.sim, method = "L-BFGS-B",lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
     }
     else{
        if(boxcox){
           lik.optim <- nlminb(ini, intgr.lik.sim.boxcox, lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
        else{
           lik.optim <- nlminb(ini, intgr.lik.sim, lower = lower.optim, upper = upper.optim,
                 fp = fixed.values, ip = ip, temp.list = temp.list)
        }
     }
  }
  ##
  if(messages.screen) 
    cat("pois.likfit.sim: end of numerical maximisation.\n")
  par.est <- lik.optim$par
  phi <- par.est[1]
  ##
  ## Values of the maximised likelihood
  ##
  if(is.R())
    loglik.max<-  - lik.optim$value
  else
    loglik.max <- - lik.optim$objective
  ##
  ## Assigning values for estimated parameters
  ##
  if(!fix.nugget.rel){
    nugget.rel <- par.est[2]
  }
  if(!fix.lambda){
    lambda <- par.est[length(par.est)]
  }
  ##
  if(is.R()) gc(verbose = FALSE)
  ##
  if(lik.type == "standard"){
     ## Computing estimated beta and sigmasq
     if((phi < 1e-12))
       siv <- list(diag.inverse = rep(1/(1+nugget.rel), n), lower.inverse = rep(0,n*(n-1)/2))
     else{
       siv <- varcov.spatial(coords = coords, cov.model = cov.model, kappa = kappa, nugget = nugget.rel, cov.pars = c(1, phi),
                   only.inv.lower.diag = TRUE)
     }
     if(boxcox){
        if(lambda == 0){
           result <- maxim.aux1(S = log(mcmc.obj$mu), invcov = siv, trend = as.vector(temp.list$xmat),
                 log.f.sim = temp.list$log.f.sim - apply(log(mcmc.obj$mu),2,sum)*(lambda-1) )
        }
        else{
           result <- maxim.aux1(S = (mcmc.obj$mu^lambda-1)/lambda,invcov = siv,trend = as.vector(temp.list$xmat),
                 log.f.sim = temp.list$log.f.sim - apply(log(mcmc.obj$mu),2,sum)*(lambda-1) )
        }
        results <- list(cov.model = cov.model, beta = result$beta, cov.pars=c(result$sigmasq, phi), nugget.rel = nugget.rel, 
                              kappa = kappa, aniso.pars = aniso.pars, lambda = lambda, loglik.max = loglik.max, call = call.fc)
     }
     else{
        result <- maxim.aux1(S = mcmc.obj$S,invcov = siv,trend = as.vector(temp.list$xmat),log.f.sim = temp.list$log.f.sim)
        results <- list(cov.model = cov.model, beta = result$beta, cov.pars = c(result$sigmasq, phi), nugget.rel = nugget.rel, 
                              kappa = kappa, aniso.pars = aniso.pars, lambda = mcmc.obj$lambda, loglik.max = loglik.max, call = call.fc)
     }
  }
  else{
     if(boxcox)
        results <- list(cov.model = cov.model, phi = phi, nugget.rel = nugget.rel, kappa = kappa, aniso.pars = aniso.pars, 
                              lambda = lambda, loglik.max = loglik.max, call = call.fc)
     else 
        results <- list(cov.model = cov.model, phi = phi, nugget.rel = nugget.rel, kappa = kappa, aniso.pars = aniso.pars, 
                              lambda = mcmc.obj$lambda, loglik.max = loglik.max, call = call.fc)
  }
  ##
  ##
  return(results)
}


