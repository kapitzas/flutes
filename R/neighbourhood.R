#' neihgbourhood
#'
#' Applies neighbourhood function to data
#'
#' @param lu land use data matrix or RasterStack object. Each column or layer represents a land use class.
#' @param cols Columns (land use classes) to apply neighbourhood function to
#' @param weights A list of weight matrices for each land use to which the neighbourhood function is applied (same length as cols)
#' @param mask a raster object containing the study area mask. NA values are outside of study area, non-NA values are within study area and must be of same length as the number of rows in the provided land use data matrix.
#' @param enr When \code{enr=TRUE}, calculated neihbourhood values are returned as enrichment factors.
#' @param format character string, default is `matrix`. When `format = "stack"`, neighbourhood layers are returned as raster stack.
#'
#' @return A matrix of neighbourhood values of the same dimensions and structure as \code{lu}.
#'
#' @import raster
#'
#' @export

neighbourhood <- function(lu, cols, weights, mask, ..., enr = F, suffix = NULL, format = "matrix"){

  inds <- which(!is.na(getValues(mask)))
  if(class(lu)[1] == "matrix"){
    if(length(inds) != nrow(lu)) stop("NA not synched")
  }
  if(class(lu)[1] == "RasterStack"){
    if(length(inds) != length(which(!is.na(lu[[1]][])))) stop("NA not synched")
    lu <- raster::as.data.frame(lu)[inds,]
  }


  out <- lu[,cols]
  for (i in 1:length(cols)){
    w <- weights[[cols[i]]]
    l <- mask

    l[inds] <- lu[,cols[i]]
    l <- raster::focal(l, w, na.rm = TRUE, pad = TRUE,...)
    if(enr == TRUE){ #Should ouput be converted to enrichment factors (according to Verburg et al 2004)
      l <- l/mean(lu[,cols[i]])
    }
    out[,i] <- getValues(l)[inds]
  }
  if(!is.null(suffix)){
    colnames(out) <- paste0(colnames(lu), "_", suffix)
  }

  if(format == "stack"){
    s <- stack(lapply(1:length(cols), function(i) mask))
    for(i in 1:nlayers(s)){
      s[[i]][inds] <- out[,i]
    }
    names(s) <- colnames(out)
    out <- s
  }
  out
}
