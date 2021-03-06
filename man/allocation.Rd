% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allocation.R
\name{allocation}
\alias{allocation}
\title{allocation}
\usage{
allocation(lu, sm, params, dmd, ln, constraint, pa = NULL)
}
\arguments{
\item{lu}{Matrix containing fractional land use data. Columns are classes, rows are cells. Each row sums to 1.}

\item{sm}{A matrix containing the predicted land use suitability model for the current time step.}

\item{params}{A list of model parameters \emph{resolution}, \emph{max_dev}, \emph{growth}, \emph{no_change} (see details).}

\item{dmd}{Matrix containing demand changes to be allocated. The first row is the current land use supply, the second row the land use demand for the next time step. Output of \code{fliutes::demand}.}

\item{ln}{A matrix containing neighbourhood variables of the current time step. These can be estimated by using the \code{flutes::neighbourhood} function.}

\item{constraint}{\code{TRUE} or \code{FALSE}. When true, constraints are applied (see details).}

\item{pa}{default is \code{pa = NULL} or raster layer containing protected areas (where the raster contains 0) from which changes are masked. All other areas (where the raster is 1) are allowed to change.}
}
\value{
A matrix containing predicted land use fractions (same format as input matrix)
}
\description{
Allocates changes in land use demand across the landscape
}
\details{
Model parameters that can be specified with \code{params} argument: \emph{max_dev} is the \% deviation allowed between allocated supply and prescribed demand. \emph{resolution} is the integer count to which fractions are converted for multinomial model. i.e. a fraction of 0.2 would become round(0.2 * resolution). \emph{growth} is the \% of the landscape (in terms of cells) where land use was 0 in a land use and contained a fraction of that land-use in the next time step, averaged across the observed time series. \emph{no_change} is a vector containing the indices of land use types that are not allowed to change at all. Applied constraints when \code{constraint = TRUE}, most cells on which land use in a class is 0 are masked from increases in that class to the next time step. Only a small subset of cells for which this is the case are allowed to change. The number of these cells is determined by \emph{growth}, which can be estimated from historic data, literature review or expert advice.
}
