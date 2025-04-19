#' Test for marginal effects of predictors on a scalar response
#'
#' @param x A numeric matrix of predictors
#' @param y A numeric vector of responses
#' @param B Number of bootstrap samples (default 199)
#' @param method Method for p-value calculation: "max", "sum", or "adaptive"
#' @return A list containing the p-value and computation time
#' @note This function is based on the C implementation by Zhang and Laber (2014) <doi:10.1080/01621459.2015.1106403>.
#' @examples
#' # Generate sample data
#' set.seed(47)
#' n <- 200
#' p <- 10
#' x <- matrix(rnorm(n*p), n, p)
#' y <- 0.25*x[,1] + rnorm(n)
#' # Run the test
#' marginal.test(x, y, B = 200, method = "adaptive")
#' marginal.test(x, y, B = 200, method = "max")
#' marginal.test(x, y, B = 200, method = "sum")
#' @export
marginal.test <- function(x, y, B = 199L, method = "adaptive")
{
  time0 <- proc.time()[3L]

  # Check input validity
  if (!is.matrix(x)) stop("x must be a matrix.")
  if (!is.numeric(x)) stop("x must be numeric.")
  if (!is.numeric(y)) stop("y must be numeric.")
  if (length(y) != nrow(x)) stop("x and y must have the same number of rows.")
  if (B <= 0) stop("B must be a positive integer.")
  if (!is.character(method) || !method %in% c("max", "sum", "adaptive")) {
    stop("Invalid method. Choose from 'max', 'sum', or 'adaptive'.")
  }
  if (B != round(B)) stop("B must be an integer.")
  if (B < 199L) warning("B is set to a low value. Consider increasing it for more reliable results.")
  if (B > 10000L) warning("B is set to a high value. This may take a long time to compute.")

  if (!is.matrix(x)) x <- cbind(x)
  n <- nrow(x)
  p <- ncol(x)

  alpha <- seq_len(p)
  num_alpha <- length(alpha)

  # Call the Rcpp implementation
  extreme <- detect_effect_rcpp(x, y, alpha, as.integer(B))

  time0 <- as.double(proc.time()[3L] - time0)

  if(method == "max") {
    return(list(
      "p_value" = extreme[1L],
      "time" = time0
    ))
  } else if(method == "sum") {
    return(list(
      "p_value" = extreme[num_alpha],
      "time" = time0
    ))
  } else if(method == "adaptive") {
    return(list(
      "p_value" = extreme[num_alpha + 1L],
      "time" = time0
    ))
  } else {
    stop("Invalid method. Choose from 'max', 'sum', or 'adaptive'.")
  }
}

