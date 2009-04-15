\name{glsm.krige}

\alias{glsm.krige}


\title{Prediction for a Generalised Linear
  Spatial Model}

\description{This function makes prediction for a
  generalised linear spatial model, using an output object from \code{\link{glsm.mcmc}}
}

\usage{
glsm.krige(mcmc.output, locations, borders, trend.l,
micro.scale=NULL, dist.epsilon= 1e-10, output)
}

\arguments{
  \item{mcmc.output}{an output file from the function
    \code{\link{glsm.mcmc}}.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data frame, or a
    list with the 2-D
    coordinates of the \eqn{N} prediction locations.  }
  \item{borders}{optional. If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this
    polygon.  }
  \item{trend.l}{specifies the trend (covariate) values at prediction
    locations. It must be of the same type as for \code{trend.d}.  }
  \item{micro.scale}{micro-scale variance. If specified, the
    nugget is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}.
    This has effect on prediction, since the the target for
    prediction is inverse link function of the ``signal'' part of \eqn{S}
    (without the measurement error part of the nugget). The
    default is \code{micro.scale = nugget}.  }
  \item{dist.epsilon}{a numeric value. Locations which are separated
    by a distance less than this value
    are considered co-located.  }
  \item{output}{parameters for controlling the output. It can take an output from \code{output.glm.control} or
    a list with elements as for the arguments in \code{output.glm.control}.
    See documentation for \code{output.glm.control}.  } 
}

\details{
  This function makes prediction for fixed parameters using an output object from
  \code{\link{glsm.mcmc}} containing the model specification and
  simulations from the posterior values of \eqn{S}{S}.
  
  The prediction consist of performing trans-Gaussian kriging on each of the simulated 
  \eqn{g^{-1}(S)}{g^{-1}(S)}-``datasets'' from the conditional
  distribution. Afterwards the predictor is obtained by taking the mean of 
  prediction means, and the prediction variance
  is obtained by taking the mean of the prediction variances plus the variance of the prediction means.
  The trans-Gaussian kriging is done by calling the function \code{krige.conv.extnd}, which is an extension of 
  \code{\link[geoR]{krige.conv}} allowing for more than one ``data
  set'', and using a second order Taylor approximation of the inverse
  link function \eqn{g^{-1}}{g^{-1}}.
}

\value{
  A list with the following components:
  
  \item{predict}{a vector with predicted values.  }
  \item{krige.var}{a vector with predicted variances.  }
  \item{mcmc.error}{estimated Monte Carlo errors on the predicted values.  }
  \item{simulations}{an \eqn{ni \times n.sim}{ni x n.sim} matrix where \eqn{ni} is the number of prediction locations and \eqn{n.sim}
    is the number of MCMC simulations. Each column 
    corresponds to a conditional simulation of the predictive
    distribution \eqn{g^{-1}(S^{*})}. Only returned if \code{sim.predict = TRUE}.  }
  \item{message}{messages about the type of prediction performed.  }
  \item{call}{the function call.  }
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://gbi.agrsci.dk/~ofch/geoRglm}.
}

\author{
  Ole F. Christensen \email{OleF.Christensen@agrsci.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{glsm.mcmc}} for MCMC simulation in a generalised linear spatial model.
}

\examples{
\dontshow{library(geoR)}
if(!exists(".Random.seed", envir=.GlobalEnv, inherits = FALSE)) set.seed(1234)
data(b50)
mcmc.5 <- mcmc.control(S.scale = 0.6, thin=1)
model.5 <- list(cov.pars=c(0.6, 0.01), beta=1, family="binomial")
outmcmc.5 <- glsm.mcmc(b50, model= model.5, mcmc.input = mcmc.5)
test2 <- glsm.krige(outmcmc.5, locations=matrix(c(0.15,0.15,0.005,0.05),2,2))
image(test2)
test3 <- glsm.krige(outmcmc.5, locations=matrix(c(0.15,0.15,0.005,0.05),2,2),
                     output=output.glm.control(sim.predict=TRUE, quantile=FALSE))
}

\keyword{spatial}
