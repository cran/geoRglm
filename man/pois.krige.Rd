\name{pois.krige}

\alias{pois.krige}
\alias{krige.conv.extnd}
\alias{mcmc.aux}
\alias{mcmc.boxcox.aux}
\alias{mcmc.pois.log}
\alias{mcmc.pois.boxcox}



\title{Conditional Simulation and Prediction for the Poisson
  Spatial model with a link function from the Box-Cox class}
\description{
  This function performs conditional simulation (by MCMC) and spatial prediction in the Poisson normal
  model with link function from the Box-Cox class for fixed covariance parameters.
  Available types of prediction are:
  \emph{sk} (simple kriging; fixed beta),
  \emph{ok} (ordinary kriging; uniform prior on beta). 
}

\usage{
pois.krige(geodata, coords = geodata$coords, data = geodata$data,
               units.m = "default", locations = NULL, borders = NULL,
               mcmc.input, krige, output)
}

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \bold{geoR} data set.
    If not provided the arguments \code{coords} and \code{data} must be provided instead.  
    The list may also contain an argument \code{units.m} as described
    below.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  } 
  \item{units.m}{ \eqn{n}-dimensional vector of observation times for the data. By default (\code{units.m = "default"}), 
    it takes \code{geodata$units.m} in case this exist and else a vector of 1's.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data frame with the 2-D
    coordinates of the \eqn{N} prediction locations.  }
  \item{borders}{optional. If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this
    polygon.  }
  \item{mcmc.input}{input parameter for the MCMC algorithm. It can take an output from \code{mcmc.control} or
    a list with elements as for the arguments in
    \code{mcmc.control}. See documentation for
    \code{mcmc.control}.\cr ATTENTION: the argument \code{S.scale} is necessary while
    all the others have default values.   }
  \item{krige}{defines the model components and the type of
    kriging. It can take an output from \code{krige.glm.control} or
    a list with elements as for the arguments in
    \code{krige.glm.control}. See documentation for \code{krige.glm.control}. } 
  \item{output}{parameters for controlling the output. It can take an output from \code{output.glm.control} or
    a list with elements as for the arguments in \code{output.glm.control}.
    See documentation for \code{output.glm.control}.  } 
}


\details{
  For simulating the conditional distribution of \eqn{S} given \eqn{y}, the Langevin-Hastings algorithm 
  with the parametrisation in Papaspilliopoulus, Roberts and Skold (2003)
is used. This algorithm is a Metropolis-Hastings algorithm, where the
proposal distribution uses gradient information from the
log-posterior distribution. 

The proposal variance (called \code{S.scale}; see \code{mcmc.control})
for the algorithm needs to be scaled
such that approximately 60 percent of the proposals are accepted. We
also recommend that the user to studies plots of the autocorrelations.
 
The prediction part of the program consist of performing trans-Gaussian kriging on each of the simulated 
\eqn{g^{-1}(S)}{g^{-1}(S)}-``datasets'' from the conditional distribution. Afterwards the predictor is obtained by taking the mean of 
prediction means, and the prediction variance
is obtained by taking the mean of the prediction variances plus the variance of the prediction means.
The trans-Gaussian kriging is done by calling the function \code{krige.conv.extnd}, which is an extension of 
\code{\link[geoR]{krige.conv}} allowing for more than one ``data
set'', and using a second order Taylor approximation of the inverse
Box-Cox transformation function \eqn{g^{-1}}{g^{-1}} when the
transformation parameter \eqn{lambda>0}; for the exponential function,
i.e. logarithmic link and \eqn{lambda=0}, an exact formula is used
instead of the Taylor approximation. 
}

\value{
  A list with the following components:
  
  \item{predict}{a vector with predicted values.  }
  \item{krige.var}{a vector with predicted variances.  }
  \item{mcmc.error}{estimated Monte Carlo errors on the predicted values.  }
  \item{beta.est}{estimate of the \eqn{\beta}{beta} parameter. Not included in the output if
    \code{type.krige = "sk"}.  }
  \item{intensity}{an \eqn{n \times n.sim}{n x n.sim} matrix with
    \eqn{n.sim} being the number of MCMC simulations. containing \eqn{g^{-1}(S_i)}. Each column 
    corresponds to a conditional simulation of the conditional distribution of \eqn{g^{-1}(S_i)} at the data locations. 
    Only returned when no prediction locations are given.  } 
  \item{acc.rate}{a data frame with acceptance rates from MCMC.
    Only returned when no prediction locations are given.  }
  \item{simulations}{an \eqn{ni \times n.sim}{ni x n.sim} matrix where \eqn{ni} is the number of prediction locations and \eqn{n.sim}
  is the number of MCMC simulations. Each column 
    corresponds to a conditional simulation of the predictive
    distribution \eqn{g^{-1}(S^{*})}. Only returned if \code{sim.predict = TRUE}.  }
  \item{call}{the function call.  }
}

\references{
  O. Papaspiliopoulus and G. O. Roberts and M. Skold
  (2003). Non-centered parameterizations for hierarchical models and
  data augmentation. \emph{Bayesian statistics 7} (eds. J. M. Bernardo,
  S. Bayarri, J. O. Berger, A. P. Dawid, D. Heckerman, A. F. M. Smith
  and M. West), Oxford University Press, 307-326.
}

Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}} for Bayesian prediction in the
  Poisson-normal model, \code{\link{binom.krige}} for prediction with
  fixed parameters in the binomial-logit normal model,
  \code{\link[geoR]{krige.conv}} for
  prediction in the linear Gaussian model.
}

\examples{
\dontshow{library(geoR)}
if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) set.seed(1234)
data(p50)
# First we scale the algorithm, and study how well the chain is mixing.
test <- pois.krige(p50, krige = list(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.2, thin = 1))
plot(log(test$intensity[45,]), type = "l")
if(!("package:stats" \%in\% search())) require("ts")
acf(log(test$intensity[45,]), type = "correlation", plot = TRUE)
\dontrun{# Now we make prediction (we decide to thin to every 10, which is the default),
# where we now use S.scale = 0.55.
test2 <- pois.krige(p50, locations = cbind(c(0.5,0.5), c(1,0.4)), 
      krige = krige.glm.control(cov.pars = c(1,1), beta = 1), 
      mcmc.input = mcmc.control(S.scale = 0.55)) }
\dontshow{
y9 <- grf(grid = expand.grid(x = seq(1, 3, l = 3), y = seq(1, 3, l = 3)), cov.pars = c(0.1, 0.2))
y9$data <- rpois(9, lambda = exp(y9$data))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1)
test2 <- pois.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = mcmc.control(S.scale = 0.5, thin = 1, n.iter=10), output = list(sim.predict = TRUE))
model2.u <- krige.glm.control(cov.pars = c(1,1), type.krige = "ok")
test2.unif.beta <- pois.krige(y9, krige = model2.u, mcmc.input =list(S.scale = 0.5, thin = 1, n.iter=10))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1, aniso.pars=c(1,2))
test2 <- pois.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = list(S.scale = 0.5, thin = 1, n.iter=10)) }
}

\keyword{spatial}













