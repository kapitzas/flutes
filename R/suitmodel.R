#' suitmodel
#'
#' Estimates parameters of multinomial, multi-class suitability model.
#'
#' @param form Model formula (only the right hand side).
#' @param lu Matrix containing land use data. Columns are classes, rows are cells. Each row sums to 1.
#' @param sub Size of subset (default `sub = NULL`) of data on which model is built.
#' @param resolution The resolution of the integer representation of proivided land use fractions (necessary in multinomial models)
#' @param data Model covariates. Each column is a covariate, rows correspond to same cells as in `lu`.
#' @param ... Further parameters passed to nnet::multinom
#'
#' @return A nnet multinomial model object
#' @details The function first converts the input matrix to integer representation by sampling from a multinomial distribution and using provided land use fractions as probability vectors. Parameters are estimated on these integers.
#' @examples
#'
#' @export



suitmodel <- function(form, lu, data, sub = NULL, resolution, decay, ...){

  #Sample random subset of data to build model on
  if(!is.null(sub)){
    subs <- sample(1:nrow(lu), sub)
  }else{
    subs <- 1:nrow(lu)
  }

  #Turn fractions into frequency counts for multinomial regression
  counts <- integerify(x = lu[subs,], resolution = resolution)
  data_sub <- as.data.frame(data[subs,])

  #Pass pre-determined formula
  f <- as.formula(paste("counts", "~", form))
  suit_model <- nnet::multinom(f, data = data_sub, decay = decay, ...)
  suit_model
}
