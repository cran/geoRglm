\name{covariog}

\alias{covariog}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Empirical Covariogram for a Model with log-link and an Underlying Gaussian Field}

\description{
  Computes the sample empirical (sample) covariogram described in Christensen, Moller and Waagepetersen (2000). 
  Output is returned as a \code{binned covariogram}.
}

\usage{
covariog(geodata, coords = geodata$coords, data = geodata$data,
         units.m = "default", uvec = "default", bins.lim = "default",
         estimator.type = c("poisson", "not-poisson"),
         max.dist = NULL, pairs.min = 2)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing elements \code{data} and \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data set.
    If not provided the arguments \code{data} and
    \code{coords} must be provided instead.
    The list may also contain an argument \code{units.m} as described below.
  } 
  \item{data}{a vector or matrix with data values.
    If a matrix is provided, each column is regarded as one variable or realization.
    Default is \code{geodata$data}, if provided.
  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the \eqn{n} data locations in each row.
    Default is \code{geodata$coords}, if provided.
  }
  \item{units.m}{ \eqn{n}-dimensional vector of observation times for the data. By default (\code{units.m = "default"}),
  it takes \code{geodata$units.m} in case this exist and else a vector of 1's.  }
  \item{uvec}{a vector with values defining the covariogram binning. The
    values of \code{uvec} defines the midpoints of the bins.\cr 
    If \eqn{uvec[1] > 0} the first bin is: \eqn{0 < u <= uvec[2] - 0.5*(uvec[2] - uvec[1])}.\cr 
    If \eqn{uvec[1] = 0} first bin is: \eqn{0 < u <= 0.5*uvec[2]},
    and \eqn{uvec[1]} is replaced by the midpoint of this interval. \cr 
    The default (\code{uvec = "default"}) is that 
    \eqn{uvec[i]=max.dist*(i-1)/14} for \eqn{i=1,\ldots,15}{i=1,...,15}. 
  }
  \item{bins.lim}{separating values for the binning. By default these values are defined via the argument of 
  \code{uvec}.  }
  \item{estimator.type}{\code{"poisson"} estimates the value \eqn{\hat{C}(0)} using the
    Poisson assumption. \code{"not-poisson"} doesn't compute \eqn{\hat{C}(0)}.
  }
  \item{max.dist}{a number defining the maximal distance for the covariogram. Pairs of locations
    separated by a larger distance than this value are ignored in the
    covariogram calculation. Default is the maximum distance
    between pairs of data locations.
  }
  \item{pairs.min}{An integer number defining the minimum number of
    pairs for the bins.
    Bins with number of pairs smaller than this
    value are ignored.
  }
}

\details{
  Covariograms can be used in geostatistical analysis for exploratory
  purposes, to estimate covariance parameters and/or to compare theoretical
  and fitted models against the empirical covariogram. 
  
  Assume that the conditional distribution of \eqn{Y_i}{Y_i} given \eqn{S_i}{S_i}
  has mean \eqn{t_i\exp(S_i)}{t_i*exp(S_i)}, where the values of \eqn{t_i} are given in \code{units.m}.
  The estimator implemented is 
  \deqn{
    \hat{C}(u) = \log\left(\frac{\frac{1}{|W_u^{\Delta}|}
      \sum_{(i,j)\in W_u^{\Delta}} Y(x_i) Y(x_j) /(t_i t_j)}{\left(\frac{1}{n}\sum_{i=1}^nY(x_i)/t_i\right)^2}\right), \ \ u > 0}{hat{C}(u) = log(frac{1/|W_u^{\Delta}|\sum_{(i,j) in W_u^{\Delta}} Y(x_i)*Y(x_j) /(t_i*t_j)}{(1/n sum_i Y(x_i)/t_i)^2}),   u > 0
  }
  When a Poisson distribution is assumed, then
  \deqn{
    \hat{C}(0) = \log\left(\frac{\frac{1}{n}\sum_{i=1}^nY(x_i)(Y(x_i)-1)/t_i^2}{\left(\frac{1}{n}\sum_{i=1}^nY(x_i)/t_i\right)^2}\right)}{
  hat{C}(0) = log(frac{1/n sum_i Y(x_i)*(Y(x_i)-1)/t_i^2}{(1/n sum_i Y(x_i)/t_i)^2}).}
}

\value{
  An object of the \code{\link[~/pack/R-1.2.1/lib/R/library/]{class}} \code{covariogram} which is a
  list with the following components:
  \item{u }{a vector with distances.  }
  \item{v }{a vector with estimated covariogram values at distances given
    in \code{u}. When \code{estimator.type = "poisson"}, the first value in \code{v} is the estimate of \eqn{\sigma^2}{sigma^2},
    \eqn{\hat{C}(0)}.  }
  \item{n }{number of pairs in each bin. When \code{estimator.type = "poisson"}, the first value in \code{n} is \code{v0}.  }
  \item{v0}{the estimate of \eqn{\sigma^2}{sigma^2}, \eqn{\hat{C}(0)}.  }
  \item{bins.lim}{Separating values for the binning provided in the function call.  }
  \item{estimator.type }{echoes the type of estimator used.  }
  \item{call}{The function call.  }  
}

\references{
  Christensen, O. F., Moller, J. and Waagepetersen R. (2000). Analysis of spatial data using generalized linear mixed models 
and Langevin-type Markov chain Monte Carlo.
\emph{Research report R-00-2009}, Aalborg University.

  Further information about \bold{geoRglm} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~christen/geoRglm}.
}

\author{
  Ole F. Christensen \email{o.christensen@lancaster.ac.uk}, \cr
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}.
}

\seealso{\code{\link{covariog.model.env}} for
  covariogram envelopes and \code{\link{plot.covariogram}} for graphical
  output. 
}

\examples{
if(is.R()) data(y50)
covar <- covariog(y50, uvec=c(1:10))
plot(covar)
## Now excluding the bin at zero (only assuming log-link).
covar2 <- covariog(y50,uvec=c(1:10), estimator.type="no")
plot(covar2)
}

\keyword{spatial}
