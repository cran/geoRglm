\name{asympvar}

\alias{asympvar}


\title{Asymptotic Variance}

\description{
Calculates the initial monotone/positive sequence estimate of the asymptotic variance from CLT (Geyer 92).
Useful for estimation of the variance of a Markov Chain Monte Carlo estimate.
}

\usage{
asympvar(timeseries, type="mon", lag.max = 100, messages)
}


\arguments{
  \item{timeseries}{a vector with a timeseries, or a matrix where the rows are different timeseries.}
  \item{type}{\code{"pos"} and \code{"mon"} gives the monotone and the positive sequence estimator, respectively, and
  \code{"all"} gives both. Default is \code{type="mon"}.  }
\item{lag.max}{maximum lag at which to calculate the asymptotic
  variance.  Default is \code{lag.max = 100}.  }
\item{messages}{logical. If \code{TRUE}, the default, status messages
    are printed while the function is running.  }
}


\value{
 A number (or a vector) with the estimate, when \code{type="mon"} or \code{type="pos"}. A list with components
 \code{mon} and \code{pos} when  \code{type="all"}
}

\references{ 
  Geyer, C. (1992). Practical Monte Carlo (with discussion). 
  \emph{Statist. Sci.} \bold{7}, 473-511.
  
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.daimi.au.dk/~olefc/geoRglm}.
}

\author{
  The orginal Splus version of this function was written by Rasmus
  Waagepetersen. R port by \cr
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}


\examples{
\dontshow{library(geoR)}
data(p50)
\dontrun{
test <- pois.krige(p50, krige = krige.glm.control(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.5, n.iter = 1000, thin = 1))
asympvar(test$intensity[45,])
ass <- asympvar(test$intensity[1:10,], type = "pos")
}
\dontshow{
test <- pois.krige(p50, krige = krige.glm.control(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.5, n.iter = 100, thin = 1))
asympvar(test$intensity[45,])
ass <- asympvar(test$intensity[1:10,], type = "pos")
}
}

\keyword{spatial}
 
 



