% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_climate.R
\name{prepare_climate}
\alias{prepare_climate}
\title{Subsets or replicate a climate data}
\usage{
prepare_climate(climate, from = "2000-04", to = "2010-11")
}
\arguments{
\item{climate}{table containing the information about monthly values for climatic data. If the climate table have exactly 12 rows it will be replicated for the number of years and months specified by \code{from} - \code{to}. Otherwise, it will be subsetted to the selected time period. If this is required, \code{year} and \code{month} columns must be included in the climate table. The minimum required columns are listed below, but additionally you can include: tmp_ave, c02, d13catm. Please refer to \code{\link{d_climate}} for example.
\itemize{
\item year: year of observation (only required for subsetting) (numeric).
\item month: months of observation (only required for subsetting) (numeric).
\item tmp_min: monthly mean daily minimum temperature (C).
\item tmp_max: monthly mean daily maximum temperature (C).
\item tmp_ave: monthly mean daily average temperature (C) (optional).
\item prcp: monthly rainfall (mm month-1).
\item rh: monthly average relative humidity (%).
\item srad: monthly mean daily solar radiation (MJ m-2 d-1).
\item frost_days: frost days per month (d month-1).
\item co2: monthly mean atmospheric co2 (ppm), required if calculate_d13c=1 (optional).
\item d13catm: monthly mean isotopic composition of air (‰), required if calculate_d13c=1 (optional).
\item n_depo: monthly nitrogen deposition (optional).
\item weibull_a: scale parameter for local wind conditions (optional).
\item weibull_k: shape parameter for local wind conditions (optional).
}}

\item{from}{year and month indicating the start of simulation. Provided in form of year-month. E.g. "2000-01".}

\item{to}{year and month indicating the end of simulation. Provided in form of year-month. E.g. "2009-12", will include December 2009 as last simulation month.}
}
\value{
a data.frame with number of rows corresponding to number of simulated month and 10 columns
}
\description{
Prepares the climate table, by either replicating the average climate for the required number of years, or by subsetting from a longer time-series of climate data.
}
\details{
This function prepares the climate table for \code{\link{run_3PG}}.

In case a user provides only average climate, this is replicated for the desired simulation period.

In case a larger climate file is provided, the simulation period is selected from this.
}
\examples{
# subsetting climate data
prepare_climate( climate = d_climate, from = '2003-04', to = '2010-11')

# replicating climate data
climate = matrix(rnorm(60), ncol = 5)
colnames(climate) = c("tmp_min", "tmp_max", "prcp", "srad", "frost_days")

prepare_climate( climate = climate, from = '2000-04', to = '2010-11')
}
\seealso{
\code{\link{run_3PG}}, \code{\link{prepare_input}}, \code{\link{prepare_parameters}}, \code{\link{prepare_sizeDist}}, \code{\link{prepare_thinning}}
}
