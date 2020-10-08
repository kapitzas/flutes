#' correlations
#'
#' Automatically reduces a data set based on estimated correlations.
#'
#' @param covs A data matrix
#' @param thresh Spearman's Rank correlation coefficient
#' @param subs number of rows to randomly sample.
#' @param mask a raster object containing the study area mask. NA values are outside of study area, non-NA values are within study area and must be of same length as the number of rows in the provided land use data matrix.
#' @param enr When \code{enr=TRUE}, calculated neihbourhood values are returned as enrichment factors (see ref).
#' 
#' @return A correlation matrix containing only selected predictors (all values hsould be < 0.8).
#' 
#' @details The function uses the Spearman's Rank Correlation coefficient. The covariate chosen to be removed from a correlated pair is the one with the higher maximum correlation coefficient with any other included covariate to maximise the amount of independent information in the final predictor set.
#'
#' @examples
#'
#' @export

correlations <- function(covs, thresh = 0.7, sub = inds) {
  subs_cor <- sample(1:nrow(covs), size = sub) #sample subset to do corr analysis on
  cors <- cor(covs[subs_cor, ], method = "spearman") #get pred correlations based on that subset
  while (min(abs(cors[abs(cors) >= thresh])) != 1){ #reduce predictor set until only uncorrelated pairs remain
    values <- cors[which(abs(cors) > thresh)]
    corellated <- which(abs(cors) > thresh)
    values[values ==1] <- NA
    corellated[which(values== max(values, na.rm = T))]
    rows_highest_cor <- which(cors == max(values, na.rm = T), arr.ind = T)[,1]
    cors_cur <- abs(cors[rows_highest_cor,])
    '%ni%' <- Negate('%in%')
    m1 <- max(cors_cur[1,][cors_cur[1,]%ni%c(max(values, na.rm = T),1)])
    m2 <- max(cors_cur[2,][cors_cur[2,]%ni%c(max(values, na.rm = T),1)])
    out <- ifelse(m1 > m2, 1, 2)
    cors <- cors[-which(colnames(cors) == names(rows_highest_cor)[out]), -which(colnames(cors) == names(rows_highest_cor)[out])]
    nrow(cors)
  }
  return(cors)
}
