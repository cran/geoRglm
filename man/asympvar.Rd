\name{asympvar}

\alias{asympvar}


\title{Asymptotic Variance}

\description{
Calculates the initial monotone/positive sequence estimate of the asymptotic variance from CLT (Geyer 92).
Useful for estimation of the variance of a Markov Chain Monte Carlo estimate.
}

\usage{
asympvar(timeseries, type="mon", lag.max = 100)
}


\arguments{
  \item{timeseries}{a vector with a timeseries, or a matrix where the rows are different timeseries.}
  \item{type}{\code{"pos"} and \code{"mon"} gives the monotone and the positive sequence estimator, respectively, and
  \code{"all"} gives both. Default is \code{type="mon"}.  }
  \item{lag.max}{maximum lag at which to calculate the asymptotic variance.  Default is \code{lag.max = 100}.  }
}


\value{
 A number (or a vector) with the estimate, when \code{type="mon"} or \code{type="pos"}. A list with components
 \code{mon} and \code{pos} when  \code{type="all"}
}

\references{ 
  Geyer, C. (1992). Practical Monte Carlo (with discussion). 
  \emph{Statist. Sci.} \bold{7}, 473-511.
  
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}


\examples{
if(is.R()) data(y50)
\testonly{library(geoR)}
test <- pois.log.krige(y50, krige = krige.glm.control(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.09, n.iter = 1000, thin = 1))
asympvar(test$intensity[45,])
ass <- asympvar(test$intensity[1:10,], type = "pos")
}

\keyword{spatial}
 
 


