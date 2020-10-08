#' suitability
#'
#' Estimates parameters of multinomial, multi-class suitability model
#'
#' @param form Model formula (only the right hand side).
#' @param lu Matrix containing land use data
#' @param sub Size of subset (default `sub = NULL`) of data on which model is built.
#' @param resolution The resolution of the integer representation of proivided land use fractions (necessary in multinomial models)
#' @param ... Further parameters passed to nnet::multinom
#'
#' @return A nnet multinomial model object
#'
#' @examples
#'
#' @export



suitmodel <- function(form, lu, data, sub = NULL, resolution,...){

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
  suit_model <- nnet::multinom(f, data = data_sub, ...)
  suit_model
}
