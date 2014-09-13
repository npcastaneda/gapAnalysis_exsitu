\name{spASD15}
\alias{spASD15}
\title{Species-specific calculation of area above 0.15 s.d.}
\description{
For a given species, loads the k-fold standard deviation rasters (thresholded and unthresholded) and computes the percentage of area above a value of 0.15
}
\usage{
spASD15(idir, spID)
}
\arguments{
  \item{idir}{directory where gap process outputs are being stored}
  \item{spID}{identifier for species}
}
\value{
A \code{data.frame} with taxon, sizeComplete, sizeCompleteUncertain, rateComplete, sizeThresholded, sizeThresholdedUncertain, rateThresholded, for the taxon in question
}
\references{
Ramirez-Villegas, J. et al. (2010) A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PLoS ONE 5(10): e13497, doi: 10.1371/journal.pone.0013497
}
\author{
The original function was developed by Julian Ramirez-Villegas for the 2010 study (Ramirez-Villegas, 2010, see reference list), and further modified by Nora Castaneda and Harold Achicanoy.
}
\seealso{
See Also as \code{\link{calcASD15}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
## The function is currently defined as
}
\keyword{ ~uncertainty }
