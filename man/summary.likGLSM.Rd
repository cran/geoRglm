\name{summary.likGLSM}

\alias{summary.likGLSM}
\alias{print.summary.likGLSM}
\alias{print.likGLSM}


\title{Summarizes Parameter Estimation Results for Generalised linear
  Spatial Models}

\description{
  Summarizes results returned by the function \code{\link{likfit.glsm}}.

  Functions are \emph{methods} for \code{\link[base]{summary}} and
  \code{\link[base]{print}} for class \code{likGLSM} and \code{summary.likGLSM}.
}

\usage{
\method{summary}{likGLSM}(object,\dots)
\method{print}{likGLSM}(x, digits = max(3, getOption("digits") - 3),\dots)
\method{print}{summary.likGLSM}(x, digits = max(3, getOption("digits") - 3),\dots)
}

\arguments{
  \item{object}{an object of \code{\link[base]{class}} \emph{likGLSM}, typically
    a result of a call to \code{\link{likfit.glsm}}.  }
  \item{x}{an object of \code{\link[base]{class}} \emph{likGLSM} or
    \code{\link[geoR]{class}} \emph{summary.likGLSM}, typically
    resulting from a call to \code{\link{likfit.glsm}}.  }
  \item{digits}{the number of significant digits to use when printing.  }
  \item{\dots}{extra arguments for \code{\link[base]{print}}.  }
}

\details{
  A detailed summary of a object of the class \code{likGLSM} is produced by
  by \code{summary.likGLSM} and printed by \code{print.summary.likGLSM}.
  This includes model specification with values of fixed and estimated parameters. 
  A simplified summary of the parameter estimation is printed by
  \code{print.likGLSM}. 
}

\value{
  \code{print.likGLSM} prints the parameter estimates and the value of the
  maximized likelihood.\cr
  \code{summary.likGLSM} returns a list with main results of a call to
  \code{\link{likfit.glsm}}.\cr
  \code{print.summary.likGLSM} prints these results on the screen (or other
  output device) in a "nice" way.
}

\references{
  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{olefc@birc.au.dk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link[geoR]{likfit}}, \code{\link[base]{print}},
  \code{\link[base]{summary}}.
}

\examples{
## See examples for the function likfit.glsm
}

\keyword{spatial}
\keyword{print}

