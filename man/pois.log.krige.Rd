\name{pois.log.krige}

\alias{pois.log.krige}
\alias{krige.glm.control}
\alias{krige.conv.extnd}
\alias{mcmc.aux}
\alias{mcmc.boxcox.aux}
\alias{mcmc.pois.log}
\alias{mcmc.pois.boxcox}



\title{Conditional Simulation and Prediction for the Poisson-log
  Spatial model}
\description{
  This function performs conditional simulation (by MCMC) and spatial prediction in the Poisson log-normal
  model for fixed covariance parameters.
  Available types of prediction are:
  \emph{SK} (simple kriging; fixed beta),
  \emph{OK} (ordinary kriging; uniform prior on beta). 
}

\usage{
pois.log.krige(geodata, coords = geodata$coords, data = geodata$data,
               units.m = "default", locations = NULL,
               mcmc.input, krige, output)
krige.glm.control(type.krige = "ok", trend.d = "cte", trend.l = "cte",
                  obj.model = NULL, beta, cov.model, cov.pars, kappa,
                  nugget, micro.scale = 0, dist.epsilon = 1e-10, 
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
  \item{krige}{defines the model components and the type of
    kriging. It can take an output from \code{krige.glm.control} or
    a list with elements as for the arguments in \code{krige.glm.control}. See section DETAILS below.  } 
  \item{mcmc.input}{input parameter for the MCMC algorithm. It can take an output from \code{mcmc.control} or
    a list with elements as for the arguments in
    \code{mcmc.control}. See documentation for
    \code{mcmc.control}.\cr ATTENTION: the argument \code{S.scale} is necessary while
    all the others have default values.   }
  \item{output}{parameters for controlling the output. It can take an output from \code{output.glm.control} or
    a list with elements as for the arguments in \code{output.glm.control}.
    See documentation for \code{output.glm.control}.  } 
  \item{type.krige}{type of prediction to be performed (minimal mean
    square error prediction). Options are
    \code{"SK"} and \code{"OK"} corresponding to prediction with fixed
    parameters (\code{type.krige = "SK"}) or prediction with a uniform
    prior on \eqn{\beta}{beta} (\code{type.krige = "OK"}). 
    Prediction using a model with covariates can be done setting \code{type.krige = "OK"} 
    and specifying the
    covariate model using the arguments \code{trend.d} and
    \code{trend.l}.  }
  \item{trend.d}{specifies the trend (covariate) values at the data
    locations.
    See documentation of \code{\link[geoR]{trend.spatial}} for
    further details.
    Default is \code{tred.d = "cte"}.   }
  \item{trend.l}{specifies the trend (covariate) values at prediction
    locations. It must be of the same type as for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{obj.model}{a list with the model parameters. Could be an
    output from \code{\link{pois.likfit.sim}}.  }
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
  \item{micro.scale}{micro-scale variance. If different from zero, the
    nugget variance is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}.  This option is not implemented.  }
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
    returned is the same scale as for the original data. Only
    \code{lambda = 0} is implemented.  }
}


\details{
  For simulating the conditional distribution of \eqn{S} given \eqn{y}, the truncated Langevin-Hastings algorithm 
(see Christensen, Moller and Waagepetersen, 2001) is used. This algorithm is a Metropolis-Hastings algorithm, where the
proposal distribution uses gradient information from the log-posterior distribution. The algorithm is described below for the case 
where \eqn{\beta}{beta} is fixed.

We use the following reparametrisation of \eqn{S}.
Let \eqn{\sigma^2V} denote the
covariance matrix of \eqn{S} and let \eqn{Q} be a Cholesky factorisation. so that 
 \eqn{V=QQ^T}.
 We can then assume that 
\eqn{ S = D\beta + \sigma Q\Gamma }, where \eqn{\Gamma}
follows a \eqn{n}-dimensional standard multivariate Gaussian
 distribution. Posterior simulations of \eqn{S} are obtained by transforming
MCMC simulations from the conditional distribution of \eqn{\Gamma} given
\eqn{Y=y}. 

The log posterior density of
 \eqn{\Gamma} given \eqn{Y=y} is
\deqn{
\log f(\gamma|y)=const(y)-\frac{1}{2}\| \gamma \|^2 +\sum_{i=1}^n y_i s_i - \exp(s_i),
}{
log f(gamma|y)=const(y)-1/2 |gamma|^2 +sum_i y_i*s_i - exp(s_i),
}
where
\eqn{(s_1,\ldots,s_n)^T=D\beta + \sigma Q \gamma}{(s_1,...,s_n)^T = D*beta + sigma*Q*gamma}.

For the truncated Langevin-Hastings update we use a truncated form of the gradient (truncating by \eqn{H_i}) of the log target density,
\deqn{
\nabla(\gamma )^{trunc}= - \gamma +\sigma Q^T\left\{y_i-\exp(s_i)\wedge H_i \right\}_{i=1}^n .
}{
nabla(gamma)^{trunc}= - gamma + sigma*Q^T*{y_i - min{exp(s_i),H_i} }_{i=1}^n .
}

The proposal \eqn{\gamma'}{gamma'} follows a multivariate normal distribution with mean vector
\eqn{\xi(\gamma)=\gamma+(h/2)\nabla(\gamma)^{trunc}}{xi(gamma)=gamma+(h/2)*nabla(gamma)^{trunc}} and covariance matrix \eqn{hI},
where \eqn{h>0} is a user-specified ``proposal variance'' (called \code{S.scale}; see \code{mcmc.control}).
  
The prediction part of the program consist of performing log-Gaussian kriging on each of the simulated 
\eqn{\exp(S)}{exp(S)}-``datasets'' from the conditional distribution. Afterwards the predictor is obtained by taking the mean of 
prediction means, and the prediction variance
is obtained by taking the mean of the prediction variances plus the variance of the prediction means.
The log-Gaussian kriging is done by calling the function \code{krige.conv.extnd}, which is an extension of 
\code{\link[geoR]{krige.conv}} allowing for more than one ``data set''. 
}

\value{
  A list with the following components:
  
  \item{predict}{a vector with predicted values.  }
  \item{krige.var}{a vector with predicted variances.  }
  \item{mcmc.error}{estimated Monte Carlo errors on the predicted values.  }
  \item{beta.est}{estimate of the \eqn{\beta}{beta} parameter. Not included in the output if
    \code{type.krige = "SK"}.  }
  \item{intensity}{an \eqn{n \times n.sim}{n x n.sim} matrix where \eqn{n.sim} is the number of MCMC simulations. Each column 
    corresponds to a conditional simulation of the conditional distribution of the random effects at the data locations. 
    Only returned when no prediction locations are given.  }
  \item{message}{messages about the type of prediction performed.  }
  \item{simulations}{an \eqn{ni \times n.sim}{ni x n.sim} matrix where \eqn{ni} is the number of prediction locations and \eqn{n.sim}
  is the number of MCMC simulations. Each column 
    corresponds to a conditional simulation of the predictive distribution. Only returned if \code{sim.predict = TRUE}.  }
  \item{message}{messages about the type of prediction performed.  }
  \item{call}{the function call.  }
}

\references{
 Christensen, O. F., Moller, J. and Waagepetersen R. (2001). Geometric ergodicity of Metropolis-Hastings algorithms for 
   conditional simulation in generalised linear mixed models. \emph{Methodology and Computing in Applied Probability}. (To appear).
  
Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{pois.krige.bayes}} for Bayesian prediction in the
  Poisson-normal model and \code{\link[geoR]{krige.conv}} for
  prediction in the linear Gaussian model.
}

\examples{
if(is.R()) data(y50)
\testonly{library(geoR)}
# First we scale the algorithm, and study how well the chain is mixing.
test <- pois.log.krige(y50, krige = list(cov.pars = c(1,1), beta = 1),
      mcmc.input = mcmc.control(S.scale = 0.09, thin = 1))
plot(log(test$intensity[45,]), type = "l")
if(is.R()) require(ts)
acf(log(test$intensity[45,]), type = "correlation", plot = TRUE)
# Now we make prediction (we decide to thin to every 20),
# where we now use S.scale = 0.07.
\dontrun{test2 <- pois.log.krige(y50, locations = cbind(c(0.5,0.5), c(1,0.4)), 
      krige = krige.glm.control(cov.pars = c(1,1), beta = 1), 
      mcmc.input = mcmc.control(S.scale = 0.07, thin = 20)) }
\testonly{
y9 <- grf(grid = expand.grid(x = seq(1, 3, l = 3), y = seq(1, 3, l = 3)), cov.pars = c(0.1, 0.2))  
y9$data <- rpois(9, lambda = exp(y9$data))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1)
test2 <- pois.log.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = mcmc.control(S.scale = 0.07, thin = 1, n.iter=10), output = list(sim.predict = TRUE))
model2.u <- krige.glm.control(cov.pars = c(1,1), beta = 1, type.krige = "ok")
test2.unif.beta <- pois.log.krige(y9, krige = model2.u, mcmc.input =list(S.scale = 0.07, thin = 1, n.iter=10))
model2 <- krige.glm.control(cov.pars = c(1,1), beta = 1, aniso.pars=c(1,2))
test2 <- pois.log.krige(y9, locations = cbind(c(0.5,0.5),c(1,0.4)),
krige = model2, mcmc.input = list(S.scale = 0.07, thin = 1, n.iter=10)) }
}

\keyword{spatial}













