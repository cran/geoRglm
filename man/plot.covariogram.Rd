\name{plot.covariogram}

\alias{plot.covariogram}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Empirical Covariogram}

\description{  
  Plots sample (empirical) covariogram computed using the
  function \code{\link{covariog}}.
}

\usage{
\method{plot}{covariogram}(x, max.dist = max(x$u), ylim = "default", type = "b",
     envelope.obj = NULL, ...)
}

%- maybe also `usage' for other objects documented here.

\arguments{ 
   \item{x}{an object of the class \code{"covariogram"}, typically an
    output of the function \code{\link{covariog}}.  }
  \item{max.dist}{maximum distance for the x-axis. The default is
    the maximum distance for which the sample covariogram was computed.}
  \item{ylim}{limits for the covariogram values in the y-axis. The default
    is from the minimum to the maximum value in \code{x$v}.  }
  \item{type}{type of line for the empirical covariogram. The default is
    \code{"b"} (dots and lines).  For further details see documentation
    for \code{\link{lines}}.  }
  \item{envelope.obj}{adds a covariogram envelope computed by
    the function \code{\link{covariog.model.env}}.  }
  \item{\dots}{other arguments to be passed to the function
    \code{\link{plot}}.  }
}

\details{
  This function allows visualisation of the empirical covariogram.
  Together with \code{\link{lines.covariomodel}} it can be used
  to compare theoretical covariogram models against the
  empirical covariogram.
}

\value{
  Produces a plot with the sample covariogram on the current graphics
  device.
  No values are returned.
}

\references{
 Further information about \bold{geoRglm} can be found at:\cr
 \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
  \code{\link{covariog}} for covariogram calculations, \code{\link{lines.covariomodel}} for
  adding lines to the current plot,
  \code{\link{covariog.model.env}} for computation of covariogram envelopes,
  and \code{\link{plot}} for generic plot function.
}

\examples{
if(is.R()) data(y50)
covario <- covariog(y50, uvec = c(1:10))  # sample covariogram
plot(covario)                      
}

\keyword{spatial}
\keyword{dplot}

