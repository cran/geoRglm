"asympvar" <- 
  function(timeseries, type = "mon", lag.max = 100)
{
  if(is.R()) require(ts) 
  if(is.vector(timeseries)) n.series <- 1
  else n.series <- nrow(timeseries) 
  if(type == "mon" | type == "all" | type == "pos") {
    if(type == "mon")
      cat(paste("calculating the initial monotone sequence estimate \n"))
    if(type == "pos") 
      cat(paste("calculating the initial positive sequence estimate \n"))
    if(type == "all") 
      cat(paste("calculating the initial positive sequence estimate, and the initial monotone sequence estimate \n"))      
  }
  else stop("Must specify type as either: mon, pos or all")
  len.Gamma <- floor(lag.max/2)-1
  if(n.series == 1){
     asy.gamma <- acf(timeseries, type = "covariance", plot = FALSE, lag.max = lag.max)$acf
     asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:len.Gamma))])
     asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:len.Gamma))])
     asy.Gamma <- asy.gamma1 + asy.gamma2
     ##--------- initial monotone sequence estimate -----------#
     kmaxpos <- min(c(which(asy.Gamma<0)-1, len.Gamma))
     if(type == "all" | type =="mon"){
       kmax <- min(c(which(diff(asy.Gamma)>0),kmaxpos))
       monvarest <- 2*sum(asy.Gamma[1:kmax])-asy.gamma[1]
       if(kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
     }
     ##--------- initial positive sequence estimate -----------#
     if(type == "pos" | type == "all"){
       posvarest <- 2*sum(asy.Gamma[1:kmaxpos])-asy.gamma[1]
       if (kmaxpos == len.Gamma) warning("value of argument lag.max is not suffiently long")
     }   
  }
  else{
     if(type == "all" | type == "pos") posvarest <- rep(1,n.series)
     if(type == "all" | type == "mon") monvarest <- rep(1,n.series)
     for(i in 1:n.series){     
        asy.gamma <- acf(timeseries[i,], type = "covariance", plot = FALSE, lag.max = lag.max)$acf
        asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:len.Gamma))])
        asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:len.Gamma))])
        asy.Gamma <- asy.gamma1 + asy.gamma2
        ##--------- initial monotone sequence estimate -----------#
        kmaxpos <- min(c(which(asy.Gamma<0)-1, len.Gamma))
        if(type == "all" | type =="mon"){
          kmax <- min(c(which(diff(asy.Gamma)>0),kmaxpos))
          monvarest[i] <- 2*sum(asy.Gamma[1:kmax])-asy.gamma[1]
          if(kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
        }
        ##--------- initial positive sequence estimate -----------#
        if(type == "pos" | type == "all"){
          posvarest[i] <- 2*sum(asy.Gamma[1:kmaxpos])-asy.gamma[1]
          if (kmaxpos == len.Gamma) warning("value of argument lag.max is not suffiently long")
        }
     }
  }
  if(type == "pos") return(posvarest)
  if(type == "all") return(list(posvarest = posvarest, monvarest = monvarest))
  if(type == "mon") return(monvarest)
}

#### consider vectorising the while loops.
