#' integerify
#'
#' Converts land use fractions into integer representation, using multinomial draws.
#'
#' @param resolution The resolution at which fractions are represented as integer (see details).
#' @param no_decrease `TRUE/FALSE` vector, flagging if a land use is allowed to decrease any where in the landscape. For example, It might be desirable to prevent decreases of Urban land to reflect the high initial investment, the cost of disbanding urban land and the high socio-economic value of urban land.
#' @param x Matrix containing data set to be integerified, typically land use data in fractional representation. Columns are classes, rows are cells. Each row sums to 1.
#' @param z Matrix containing initial land use supply before iterations begin. This option is required in the `allocation` function to deduct the supply from classes that are not supposed to decrease from the total supply to be allocated in each iteration.
#' @return A matrix containing land use fractions in integer representation. Each row sums to `resolution`.
#'
#' @details Integers are sampled from a multinomial distribution using land use fractions as probability vector. Instead of summing to 1, sampled integers on each cell sum to `resolution`. For example, when `resolution = 10000` a fraction of 0.2 would approximate a value of 2000, including some very small stochasticity from the multinomial draw.
#' @export


#Integerify function
integerify <- function (x, z = NULL,  resolution = resolution, no_decrease = NULL) {

  n <- nrow(x)
  k <- ncol(x)

  if (any(no_decrease)) { #if there is at least one class that deosn't decrease (typically urban)
    min_allocation <- matrix(0, n, k)
    for (class in which(no_decrease)) { #loop through classes
      min_allocation[, class] <- z[, class]
    }
  } else {
    min_allocation <- NULL
  }

  # if a minimum number must be placed in a certain column, withold them and just allocate the others
  if (!is.null(min_allocation)) {
    to_allocate <- resolution - rowSums(min_allocation)

    # need to reduce the probability of allocating the remainder to these cells too
    # so make x sum to resolution, subtract the minima, and make sure the result is positive
    x <- x - min_allocation #might have to be other way round
    x <- resolution * x / rowSums(x) #* resolution because gets scaled to 0-1

    x <- pmax(x, 0)

  } else {
    to_allocate <- rep(resolution, n)
  }

  # random draws of the non-reserved ones
  x <- extraDistr::rmnom(n, to_allocate, prob = x)

  # add on the reserved ones again
  if(!is.null(min_allocation)) {
    x <- x + min_allocation
  }
  x
}
