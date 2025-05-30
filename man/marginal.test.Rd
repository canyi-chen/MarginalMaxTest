% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/marginal.R
\name{marginal.test}
\alias{marginal.test}
\title{Test for marginal effects of predictors on a scalar response}
\usage{
marginal.test(x, y, B = 199L, method = "adaptive")
}
\arguments{
\item{x}{A numeric matrix of predictors}

\item{y}{A numeric vector of responses}

\item{B}{Number of bootstrap samples (default 199)}

\item{method}{Method for p-value calculation: "max", "sum", or "adaptive"}
}
\value{
A list containing the p-value and computation time
}
\description{
Test for marginal effects of predictors on a scalar response
}
\note{
This function is based on the C implementation by Zhang and Laber (2014) \url{doi:10.1080/01621459.2015.1106403}.
}
\examples{
# Generate sample data
set.seed(47)
n <- 200
p <- 10
x <- matrix(rnorm(n*p), n, p)
y <- 0.25*x[,1] + rnorm(n)
# Run the test
marginal.test(x, y, B = 200, method = "adaptive")
marginal.test(x, y, B = 200, method = "max")
marginal.test(x, y, B = 200, method = "sum")
}
