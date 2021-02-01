#' neihgbourhood
#'
#' Applies neighbourhood function to data
#'
#' @param lu land use data matrix. Each column represents a land use class, each row a cell in the landscape
#' @param cols Columns (land use classes) to apply neighbourhood function to
#' @param weights A list of weight matrices for each land use to which the neighbourhood function is applied (same length as cols)
#' @param mask a raster object containing the study area mask. NA values are outside of study area, non-NA values are within study area and must be of same length as the number of rows in the provided land use data matrix.
#' @param enr When \code{enr=TRUE}, calculated neihbourhood values are returned as enrichment factors.
#'
#' @return A matrix of neighbourhood values of the same dimensions and structure as \code{lu}.
#'
#' @examples
#'
#' @export

neighbourhood <- function(lu, cols, weights, mask, ..., enr = F, suffix = NULL){
  if(length(which(!is.na(mask[]))) != nrow(lu)) stop("mask and land use not same length")
  out <- lu[,cols]
  for (i in 1:length(cols)){
    w <- weights[[cols[i]]]
    l <- mask
    inds <- which(!is.na(getValues(l)))
    l[inds] <- lu[,cols[i]]
    l <- raster::focal(l, w, na.rm = TRUE, pad = TRUE,...)
    if(enr == TRUE){ #Should ouput be converted to enrichment factors (according to Verburg et al 2004)
      l <- l/mean(lu[,cols[i]])
    }
    out[,i] <- getValues(l)[inds]
  }
  if(!is.null(suffix)){
    colnames(out) <- paste0("lu_", cols, "_", suffix)
  }
  out
}
