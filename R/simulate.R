#' Simulate a multiway contingency table
#'
#' @description Generates a contingency table with the frequencies generated
#' by taking a multinomial random sample from a randomly
#' generated parameter vector `beta`.
#'
#' @param nvars number of variables.
#' @param levels the number of levels of each of the variables.
#' @param n population size.
#' @param betas vector with values from which to draw the `beta` parameters.
#' @param seed seed for reproducibility.
#'
#' @return A lists with:
#' * `data` data frame with the observed data.
#' * `beta` vector with the true beta parameters.
#' * `fmax` formula of the maximal model.
#'
#' @examples
#' simdat(seed = 1)
#'
#' @importFrom stats setNames rmultinom
#' @importFrom dplyr %>% mutate
#' @importFrom utils combn
#' @export

simdat  <- function(nvars = 3, levels = 3, n = 1000, betas = seq(-1, 1, by = .5), seed = NULL)
{

  if (!is.null(seed)) set.seed(seed)

  # make the data excluding the frequencies
  d  <- expand.grid(lapply(1:nvars, \(x)x = factor(1:levels)), KEEP.OUT.ATTRS = F)
  d <- setNames(d, LETTERS[1:nvars])

  # rhs of formula for the saturated model
  fsat   <- formula(paste("~", paste(LETTERS[1:nvars], collapse = "*")))

  # design matrix saturated model
  X      <- model.matrix(fsat, d)
  orders <- sapply(colnames(X)[-1], \(x)nchar(x)) %>%
    rank(ties.method = "min") %>%
           factor(labels = 1:nvars) %>%
           as.numeric()

  # draw the beta parameters with smaller values for higher interactions
  beta <- sample(betas, ncol(X) - 1, replace = T) / sqrt(orders)

  # add intercept
  beta    <- c(-log(sum(exp(X[, -1] %*% beta))) + log(n), beta)
  beta    <- setNames(beta, colnames(X))

  # draw observed from expected frequencies and add to data
  m   <- exp(X %*% beta)
  Freq <- c(rmultinom(1, n, prob = m / n))
  data <- mutate(d, Freq = Freq)

  modelterms <- attr(terms(fsat, data = data), "term.labels")
  suffstats  <- crossprod(X, data$Freq)
  suff0      <- names(suffstats[suffstats == 0, 1])
  hier0      <- sapply(suff0, \(x) strsplit(x, ":"))

  effects    <- attr(terms(fsat, data = data), "term.labels")

  if (length(suff0) > 2)
  {
    for (i in length(suff0):1)
    {
      k          <- length(hier0[[i]])
      comb       <- combn(hier0[[i]], k - 1)
      subsuffs   <- apply(comb, 2, \(x)paste0(x, collapse = ":"))
      candidates <- names(hier0[lapply(hier0, length) == k - 1])

      if (any(candidates %in% subsuffs))
      {
        effect  <- gsub("[1-9]", "", suff0[[i]])
        effects <- effects[!effects %in% effect]
      }
    }
  }

  fmax <- formula(paste("Freq ~", paste(effects, collapse = "+")))

  list(beta = beta,
       data = data,
       fmax = fmax)
}



