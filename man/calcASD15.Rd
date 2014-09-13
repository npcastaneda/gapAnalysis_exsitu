\name{calcASD15}
\alias{calcASD15}
\title{All species calculation of area above 0.15 s.d.}
\description{
For all species under analysis, runs the function \code{spASD15} and returns a complete \code{data.frame} with all species information
}
\usage{
calcASD15(idir)
}
\arguments{
  \item{idir}{directory where gap process outputs are being stored}
}
\value{
A \code{data.frame} with taxon, sizeComplete, sizeCompleteUncertain, rateComplete, sizeThresholded, sizeThresholdedUncertain, rateThresholded, for all taxa being analysed. The function will also write a .csv file

}
\references{
Ramirez-Villegas, J. et al. (2010) A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PLoS ONE 5(10): e13497, doi: 10.1371/journal.pone.0013497
}
\author{
The original function was developed by Julian Ramirez-Villegas for the 2010 study (Ramirez-Villegas, 2010, see reference list), and further modified by Nora Castaneda and Harold Achicanoy.}
\seealso{
See Also as \code{\link{spASD15}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
## The function is currently defined as
}
\keyword{ ~uncertainty }
