% Generated by roxygen2 (4.0.2): do not edit by hand
\name{gibbsHMM_PT_wDist2}
\alias{gibbsHMM_PT_wDist2}
\title{gibbsHMM_PT}
\usage{
gibbsHMM_PT_wDist2(YZ, M = 2000, K = 10, alphaMAX = 1, type = 1,
  alphaMin = 0.001, J = 20, densTrue, SuppressAll = "FALSE")
}
\arguments{
\item{x,}{alpha, log=False}
}
\description{
parallel tempering with a column prior - option to mix over column or stick to j=1
}
\examples{
dDirichlet(c(.1, .9), c(0.1,0.1))
}
\keyword{dirichlet}

