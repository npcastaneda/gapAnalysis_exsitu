\name{prepareMasks}
\alias{prepareMasks}

\title{
Converts the masking shapefile into a raster and creates an additional raster for measuring areas.
}
\description{
Converts the \code{shapefile} representing the complete geographic area under study into a \code{raster}, and prepares a dummy raster that will help calculating areas in square kilometers when the Geographic Representativeness Score (GRS) is estimated.
}
\usage{
createMasks(shpMask,env_dir)
}

\arguments{
  \item{shpMask}{\code{shapefile} name that represents the complete geographic extent of all the taxa to be analyzed}
  \item{env_dir}{directory containing the bioclim \code{rasters} that will be used for preparing the species distribution models }
}
\value{
  \item{mask}{\code{raster} representing the complete geographic extent of the taxa under analysis}
  \item{cellArea}{\code{raster} that will be used to calculate the approximate area in square kilometers in a Raster object}
}
\references{
Ramirez-Villegas, J. et al. (2010) A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PLoS ONE 5(10): e13497, doi: 10.1371/journal.pone.0013497
}
\author{
The original function was developed by Julian Ramirez-Villegas for the 2010 study (Ramirez-Villegas, 2010, see reference list), and further modified by Nora Castaneda and Harold Achicanoy.
}

\seealso{
See Also as \code{\link{sizeDR}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}

