\name{splitHG}
\alias{splitHG}
\title{Splits data into germplasm accessions (G) and reference samples (H)}
\description{
Takes the main file containing occurrences, creates two separate files (G and H) and maps the density of samples
}
\usage{
splitHG(occ, crop, lon, lat)
}
\arguments{
\item{occ}{A dataframe containing the records for one or more species}

\item{crop}{character object containing the name of crop gene pool under analysis}

\item{lon}{character object. Name of the field in occ containing the longitude}

\item{lat}{character object. Name of the field in occ containing the latitude}
}

\value{
A \code{data.frame} with taxon, latitude, longitude and a field explicitly defining if each record corresponds to a germplasm accessions (G) or a reference sighting (H).
}
\references{
Ramirez-Villegas, J. et al. (2010) A gap analysis methodology for collecting crop genepools: a case study with Phaseolus beans. PLoS ONE 5(10): e13497, doi: 10.1371/journal.pone.0013497
}
\author{
The original function was developed by Julian Ramirez-Villegas for the 2010 study (Ramirez-Villegas, 2010, see reference list), and further modified by Nora Castaneda and Harold Achicanoy.
}
\seealso{
See Also as \code{\link{countRecords}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--  or do  help(data=index)  for the standard data sets.
## The function is currently defined as
}
