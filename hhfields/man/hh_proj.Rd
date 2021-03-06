% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hh_decomposition.R
\name{hh_proj}
\alias{hh_proj}
\title{Perform Hodge-Helmholtz Decomposition}
\usage{
hh_proj(dx, dy, Xi, Xb, Nb, Xe, kerngen)
}
\arguments{
\item{dx}{delta x}

\item{dy}{delta y}

\item{Xi}{(matrix: (n_interior,2)) X interior nodes}

\item{Xb}{(matrix: (n_boundary, 2)) X boundary nodes}

\item{Nb}{(matrix: (n_boundary, 2)) Normal vectors at boundary nodes}

\item{Xe}{(matrix: (n, 2)) X domain points at which to evaluate projection}

\item{kerngen}{kernel generator function}
}
\value{
List containing 
\itemize{
\item{coef}{projection coefficients}
\item{divfree}{divergence free projection at Xe}
\item{curlfree}{curlfree projection at Xe}
}
}
\description{
Perform Hodge-Helmholtz Decomposition
}
