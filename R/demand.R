#' land use demand
#'
#' Extracts aggregated (landscape-wide average)  land use changes from observed land use time series for model validation.
#'
#' @param landuse Matrix containing observed land use time series. Rows are cells, columns are land use classes + time steps.
#' @param ts Output time steps. When `ts` is longer than the length of the provided time series, demands for missing time steps are interpolated.
#' @param k Number of land use classes.
#' @param type Should the landscape mean fraction be calculated (option _mean_) or the number of cells changing in each class counted (option _sum_)?
#' @param path Default `path = NULL`. Output path to write result to (.rds format).
#' @return A matrix containing observed and interpolated average demand for different land use classes. Each column corresponds to a land use class, each row to a time step.
#'
#' @export

demand <- function(landuse, ts, path = NULL, k, type = "mean"){
  lu <- landuse

  #Get the current land use supply in each class
  if(type == "mean"){
    class_supply <- as.numeric(colMeans(lu))
  }

  if(type == "sum"){
    class_supply <- colSums(lu)
  }

  years <- ts[1]:tail(ts,1)
  demand <- matrix(NA, nrow = length(years), ncol = k + 2)
  demand[,1] <- years
  obs_ind <- which(demand[,1]%in%ts)

  for(i in 1:(ncol(lu)/k)){
    demand[obs_ind[i], 2:(k+1)] <- class_supply[seq((i * k - (k-1)), i * k, by = 1)]
  }

  #interpolate intended time steps
  for(i in 1:k){
    demand[, i + 1] <- approx(demand[,1], demand[,i + 1], xout = years)$y
  }

  demand[, k + 2] <- rowSums(demand[,-1], na.rm = TRUE)
  if(!is.null(path)){
    saveRDS(demand, path)
  }
  return(demand)
}
