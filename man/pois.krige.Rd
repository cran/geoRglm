\name{pois.krige}

\alias{pois.krige}
\alias{krige.glm.control}
\alias{krige.glm.check.aux}
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
krige.glm.control(type.krige = "sk", trend.d = "cte", trend.l = "cte",
                  obj.model = NULL, beta, cov.model, cov.pars, kappa,
                  nugget, micro.scale, dist.epsilon = 1e-10, 
                  aniso.pars, lambda)
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
    \code{krige.glm.control}.  } 
  \item{output}{parameters for controlling the output. It can take an output from \code{output.glm.control} or
    a list with elements as for the arguments in \code{output.glm.control}.
    See documentation for \code{output.glm.control}.  } 
  \item{type.krige}{type of prediction to be performed (minimal mean
    square error prediction). Options are
    \code{"sk"} and \code{"ok"} corresponding to prediction with fixed
    parameters (\code{type.krige = "sk"}), which is the default, or prediction with a uniform
    prior on \eqn{\beta}{beta} (\code{type.krige = "ok"}). 
    Prediction using a model with covariates can be done by specifying the
    covariate model using the arguments \code{trend.d} and
    \code{trend.l}.  }
  \item{trend.d}{specifies the trend (covariate) values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details.
    Default is \code{trend.d = "cte"}.   }
  \item{trend.l}{specifies the trend (covariate) values at prediction
    locations. It must be of the same type as for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{obj.model}{a list with the model parameters.  }
  \item{beta}{numerical value of the mean (vector) parameter.
    Only used if \code{type.krige="sk"}.  }
  \item{cov.model}{string indicating the name of the model for the
    correlation function. Further details  in the
    documentation for \code{\link[geoR]{cov.spatial}}.  }
  \item{cov.pars}{a vector with the 2 covariance parameters \eqn{\sigma^2}{sigma^2},
    and \eqn{\phi}{phi} for the underlying Gaussian field.  }
  \item{kappa}{additional smoothness parameter required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} and
    \code{"gneiting.matern"}.   }
  \item{nugget}{the value of the nugget parameter
    \eqn{\tau^2}{tau^2} for the underlying Gaussian field. Default is
    \code{nugget = 0}.  }
  \item{micro.scale}{micro-scale variance. If specified, the
    nugget is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}.
    This has effect on prediction where the ``signal'' part of \eqn{S}
    (without the measurement error part of the nugget) is predicted. The
    default is \code{micro.scale = nugget}.  }
  \item{dist.epsilon}{a numeric value. Locations which are separated
    by a distance less than this value
    are considered co-located.  }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link[geoR]{coords.aniso}}.  }
  \item{lambda}{numeric value of the Box-Cox transformation parameter.
    The value \eqn{\lambda = 1}{lambda = 1} corresponds to
    no transformation and \eqn{\lambda = 0}{lambda = 0} corresponds to
    the log-transformation. 
    Prediction results are back-transformed and
    returned is the same scale as for the original data.  }
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
set'', and for \eqn{lambda>0} a second order Taylor approximation. 
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
  \item{message}{messages about the type of prediction performed.  }
  \item{call}{the function call.  }
}

\references{
  O. Papaspiliopoulus and G. O. Roberts and M. Skold
  (2003). Non-centered parameterizations for hierarchical models and
  data augmentation. \emph{Bayesian statistics 7} (eds. J. M. Bernardo,
  S. Bayarri, J. O. Berger, A. P. Dawid, D. Heckerman, A. F. M. Smith
  and M. West), Oxford University Press (to appear).
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
\testonly{library(geoR)}
if(!exists(".Random.seed")) set.seed(1234)
data(p50)
# First we scale the algorithm, and study how well the chain is mixing.
test <- pois.krige(p50, krige = list(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.2, thin = 1))
plot(log(test$intensity[45,]), type = "l")
require(ts)
acf(log(test$intensity[45,]), type = "correlation", plot = TRUE)
# Now we make prediction (we decide to thin to every 10, which is the default),
# where we now use S.scale = 0.55.
\dontrun{test2 <- pois.krige(p50, locations = cbind(c(0.5,0.5), c(1,0.4)), 
      krige = krige.glm.control(cov.pars = c(1,1), beta = 1), 
      mcmc.input = mcmc.control(S.scale = 0.55)) }
\testonly{
y9 <- grf(grid = expand.grid(x = seq(1, 3, l = 3), y = seq(1, 3, l = 3)), cov.pars = c(0.1, 0.2))
y9$data <- rpois(9, lambda = exp(y9$data))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1)
test2 <- pois.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = mcmc.control(S.scale = 0.5, thin = 1, n.iter=10), output = list(sim.predict = TRUE))
model2.u <- krige.glm.control(cov.pars = c(1,1), beta = 1, type.krige = "ok")
test2.unif.beta <- pois.krige(y9, krige = model2.u, mcmc.input =list(S.scale = 0.5, thin = 1, n.iter=10))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1, aniso.pars=c(1,2))
test2 <- pois.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = list(S.scale = 0.5, thin = 1, n.iter=10)) }
}

\keyword{spatial}













