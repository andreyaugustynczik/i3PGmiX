#' @title Subsets or replicate a climate data
#' @description Prepares the climate table, by either replicating the average climate for the required number of years, or by subsetting from a longer time-series of climate data.
#'
#' @param climate  table containing the information about monthly values for climatic data. If the climate table have exactly 12 rows it will be replicated for the number of years and months specified by \code{from} - \code{to}. Otherwise, it will be subsetted to the selected time period. If this is required, \code{year} and \code{month} columns must be included in the climate table. The minimum required columns are listed below, but additionally you can include: tmp_ave, c02, d13catm. Please refer to \code{\link{d_climate}} for example.
#' \itemize{
#' \item year: year of observation (only required for subsetting) (numeric).
#' \item month: months of observation (only required for subsetting) (numeric).
#' \item tmp_min: monthly mean daily minimum temperature (C).
#' \item tmp_max: monthly mean daily maximum temperature (C).
#' \item tmp_ave: monthly mean daily average temperature (C) (optional).
#' \item prcp: monthly rainfall (mm month-1).
#' \item rh: monthly average relative humidity (%).
#' \item srad: monthly mean daily solar radiation (MJ m-2 d-1).
#' \item frost_days: frost days per month (d month-1).
#' \item co2: monthly mean atmospheric co2 (ppm), required if calculate_d13c=1 (optional).
#' \item d13catm: monthly mean isotopic composition of air (‰), required if calculate_d13c=1 (optional).
#' \item n_depo: monthly nitrogen deposition (optional).
#' \item weibull_a: scale parameter for local wind conditions (optional).
#' \item weibull_k: shape parameter for local wind conditions (optional).
#' }
#' @param from year and month indicating the start of simulation. Provided in form of year-month. E.g. "2000-01".
#' @param to  year and month indicating the end of simulation. Provided in form of year-month. E.g. "2009-12", will include December 2009 as last simulation month.
#'
#' @details This function prepares the climate table for \code{\link{run_3PG}}.
#'
#' In case a user provides only average climate, this is replicated for the desired simulation period.
#'
#' In case a larger climate file is provided, the simulation period is selected from this.
#'
#' @return a data.frame with number of rows corresponding to number of simulated month and 10 columns
#'
#' @seealso \code{\link{run_3PG}}, \code{\link{prepare_input}}, \code{\link{prepare_parameters}}, \code{\link{prepare_sizeDist}}, \code{\link{prepare_thinning}}
#'
#' @example inst/examples/prepare_climate-help.R
#'
#' @export
#'
prepare_climate <- function(
  climate,
  from = '2000-04',
  to = '2010-11'
){

  # Test for the columns consistensy
  if( !all(c("tmp_min","tmp_max","prcp","srad","frost_days") %in% colnames(climate)) ){
    stop( 'Climate table must include the following columns: tmp_min, tmp_max, prcp, srad, frost_days' )
  }

  # Test for NA
  if( any( is.na(climate[c("tmp_min","tmp_max","prcp","srad","frost_days")]) ) ){
    stop( "Climate table should not contain NA's" )
  }

  # prepare the time period
  from = as.Date(paste(from,"-01",sep=""))
  to = as.Date(paste(to,"-01",sep=""))

  if( from >= to ){
    stop( 'The start date is later than the end date' )
  }

  # make data.frame
  climate = data.frame(climate)

  # Replicate or subset the data
  if( dim(climate)[1] == 12 ){

    n_years <- as.numeric(format(to,'%Y')) - as.numeric(format(from,'%Y')) + 1
    month_i <- as.numeric(format(from,'%m'))
    month_e <- as.numeric(format(to,'%m'))

    climate = do.call("rbind", replicate(n_years, climate, simplify = FALSE))
    climate$year = rep( as.numeric(format(from,'%Y')):as.numeric(format(to,'%Y')), each = 12)
    climate$month = rep(1:12, times = n_years)

    if( month_i > 1 ){
      climate = climate[-c(1:(month_i-1)),]
    }

    if( month_e < 12 ){
      climate = climate[1:(nrow(climate)-(12-month_e)),]
    }

  } else {

    # test if year and month column are present
    if( !all(c("year", "month") %in% colnames(climate)) ){
      stop( 'Climate table must include year and month for subsettins.' )
    }

    climate$date = as.Date( paste(climate$year, '-', climate$month, "-01",sep="") )

    # Test if we climate data cover the requested range
    if( any( from < min(climate$date), to > max(climate$date)) ){
      stop( 'Requested time period is outside of providate dates in climate table.')
    }

    climate = climate[climate$date >= from & climate$date <= to, ]

  }

  # Add Average temperature if missing
  if( !'tmp_ave' %in% colnames(climate) ){
    climate$tmp_ave = (climate$tmp_min + climate$tmp_max) / 2
  }

  # Add relative humidity if missing
  if( !'rh' %in% colnames(climate) ){
    climate$rh = 0
  }

  # Add VPD if missing
  if( !'vpd_day' %in% colnames(climate) ){
    climate$vpd_day = get_vpd( climate$tmp_min, climate$tmp_max, climate$rh)
  }


  # Add default CO2 if missing
  if( !'co2' %in% colnames(climate) ){
    climate$co2 = 350
  }

  # Add default d13catm if missing
  if( !'d13catm' %in% colnames(climate) ){
    climate$d13catm = -7.1
  }

  # Add default nitrogen deposition if missing
  if( !'n_depo' %in% colnames(climate) ){
    climate$n_depo = 0
  }

  # Add default wind if missing
  if( !'weibull_a' %in% colnames(climate) ){
    climate$weibull_a = 0
  }

  # Add default wind if missing
  if( !'weibull_k' %in% colnames(climate) ){
    climate$weibull_k = 0
  }

  # Select final table
  climate = climate[,c("year", "month",'tmp_min', 'tmp_max', 'tmp_ave', 'prcp','rh', 'srad', 'frost_days', 'vpd_day', 'co2', 'd13catm','n_depo','weibull_a','weibull_k')]

  return( climate )
}




get_vpd <- function(tmin, tmax, rh){
  # internal function to calculate VPD if not available
  if (max(rh) > 0){
    vpd_min = 0.610780 * exp(17.2690 * tmin / (237.30 + tmin))
    vpd_max = 0.610780 * exp(17.2690 * tmax / (237.30 + tmax))
  # Convert kPa to mbar
    vpd_day = 10 * (0.5 * (vpd_max + vpd_min) * (1 - rh / 100))
  } else {
    vpd_min = 6.10780 * exp(17.2690 * tmin / (237.30 + tmin))
    vpd_max = 6.10780 * exp(17.2690 * tmax / (237.30 + tmax))

    vpd_day = (vpd_max - vpd_min) / 2
  }


  return(vpd_day)
}
