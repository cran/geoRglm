\name{lines.covariomodel}

\alias{lines.covariomodel}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Line with a  Covariogram Model }
\description{
  This function adds a theoretical covariogram to the
  current plot. The covariogram model to be
  added is typically with parameters estimated from the data.
}

\usage{
\method{lines}{covariomodel}(x, max.dist = x$max.dist, ...)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{covariomodel} which is a list containing
    information about the model parameters.  }
  \item{max.dist}{maximum distance (x-axis) to compute and draw the line
    representing the covariogram model.
    The default is the distance given by \code{x$max.dist}.  }
  \item{n.points}{number of points between \code{0} and \code{max.dist} to
    compute the values of the covariogram model.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{lines}}.  }  
}

\details{
  Allows theoretical covariogram(s) to be added to a plot.
  Together with \code{\link{plot.covariogram}} can be
  used to compare sample covariograms against fitted models.
}

\value{
  A line with a covariogram model is added to a plot on the current
  graphics device. No values are returned.
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{plot.covariogram}}, \code{\link{lines}}.}

\examples{
\testonly{library(geoR)}
sim <- grf(100, cov.pars = c(0.1, 0.2))
sim$data <- rpois(100, lambda = exp(sim$data+1))
# data generated from the poisson-log normal model
covario <- covariog(sim, max.dist = 1)  # sample covariogram
plot(covario)                      
parmval <- list(cov.model = "exponential", cov.pars = c(0.1, 0.1), 
                max.dist = 1)
class(parmval) <- "covariomodel"
lines(parmval, lty = 2)
}

\keyword{spatial}
\keyword{aplot}
