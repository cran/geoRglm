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
  if(n.series == 1){
     asy.gamma <- acf(timeseries, type = "covariance", plot = FALSE, lag.max = lag.max)$acf
     asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:(length(asy.gamma)/2 - 1)))])
     asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:(length(asy.gamma)/2 - 1)))])
     asy.Gamma <- asy.gamma1 + asy.gamma2
     len.Gamma <- length(asy.Gamma)
     ##--------- initial monotone sequence estimate -----------#
     varest <- asy.gamma[1] + 2 * asy.gamma[2]
     kmax <- 2
     while(asy.Gamma[kmax] > asy.Gamma[kmax + 1] & (asy.Gamma[kmax] > 0) & kmax < len.Gamma) {
       varest <- varest + 2 * asy.Gamma[kmax]
       kmax <- kmax + 1
       if (kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
     }
     if(type == "pos" | type == "all") {
       ##--------- initial positive sequence estimate -----------#
       cat(paste("calculating the initial positive sequence estimate"))
       posvarest <- asy.gamma[1] + 2 * asy.gamma[2]
       while((asy.Gamma[kmax] > 0) & kmax < len.Gamma) {
         posvarest <- posvarest + 2 * asy.Gamma[kmax]
         kmax <- kmax + 1
         if (kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
       }
       if(type == "pos")
         return(posvarest)
       else return(list(posvarest = posvarest, monvarest = varest))
     }
     else return(varest)
  }
  else{
     if(type == "all" | type == "pos") posvarest <- rep(1,n.series)
     if(type == "all" | type == "mon") monvarest <- rep(1,n.series)
     for(i in 1:n.series){     
        asy.gamma <- acf(timeseries[i,], type = "covariance", plot = FALSE, lag.max = 100)$acf
        asy.gamma1 <- c(asy.gamma[(1 + 2 * c(0:(length(asy.gamma)/2 - 1)))])
        asy.gamma2 <- c(asy.gamma[(2 + 2 * c(0:(length(asy.gamma)/2 - 1)))])
        asy.Gamma <- asy.gamma1 + asy.gamma2
        len.Gamma <- length(asy.Gamma)
        ##--------- initial monotone sequence estimate -----------#
        varest <- asy.gamma[1] + 2 * asy.gamma[2]
        kmax <- 2
        while(asy.Gamma[kmax] > asy.Gamma[kmax + 1] & (asy.Gamma[kmax] > 0) & kmax < len.Gamma) {
          varest <- varest + 2 * asy.Gamma[kmax]
          kmax <- kmax + 1
          if (kmax == len.Gamma) warning("value of argument lag.max is not suffiently long")
        }
        if(type == "pos" | type == "all") {
          posvarest[i] <- asy.gamma[1] + 2 * asy.gamma[2]
          while((asy.Gamma[kmax] > 0) & kmax < len.Gamma) {
            posvarest[i] <- posvarest[i] + 2 * asy.Gamma[kmax]
            kmax <- kmax + 1
            if (kmax == len.Gamma)warning("value of argument lag.max is not suffiently long")
          }
        }
        if(type == "all" | type =="mon") monvarest[i] <- varest
     }
  }
  if(type == "pos") return(posvarest)
  if(type == "all") return(list(posvarest = posvarest, monvarest = varest))
  if(type == "mon") return(monvarest)     
}


