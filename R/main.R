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
#' * `seed` the seed used for the data simulation.
#'
#' @examples
#' simdat(seed = 1)
#'
#' @importFrom stats setNames rmultinom
#' @importFrom dplyr %>% mutate
#' @importFrom utils combn
#' @export

simdat  <- function(nvars = 3, levels = 3, n = 1000, betas = -1:1, seed = sample(1e+5, 1))
{

  set.seed(seed)

  # make the data excluding the frequencies
  d      <- expand.grid(lapply(1:nvars, \(x)x = factor(1:levels)), KEEP.OUT.ATTRS = F)
  d      <- setNames(d, LETTERS[1:nvars])

  # rhs of formula for the saturated model
  fsat   <- formula(paste("~", paste(LETTERS[1:nvars], collapse = "*")))

  # design matrix saturated model
  X      <- model.matrix(fsat, d)
  orders <- sapply(colnames(X)[-1], \(x)nchar(x)) %>%
    rank(ties.method = "min") %>%
    factor(labels = 1:nvars) %>%
    as.numeric()

  # draw the beta parameters with smaller values for higher interactions
  beta       <- sample(betas, ncol(X) - 1, replace = T) / sqrt(orders)

  # add intercept
  beta       <- c(-log(sum(exp(X[, -1] %*% beta))) + log(n), beta)
  beta       <- setNames(beta, colnames(X))

  # draw observed from expected frequencies and add to data
  m          <- exp(X %*% beta)
  Freq       <- c(rmultinom(1, n, prob = m / n))
  data       <- mutate(d, Freq = Freq)

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
       fmax = fmax,
       seed = seed)
}

#' Model selection with step and L0 penalty based on the BIC
#'
#' @description
#' Model selection with step and L0 penalty based on the BIC
#'
#' @param object object made with \code{\link{simdat}}.
#' @param lambdarange range of the lambda values for the L0 penalties.
#' @param B number of lambda values in `lambdarange`.
#'
#' @return
#' A list with elements
#' * `out` data frame containing:
#'    * `method` model selection method
#'    * `ic` information criterion (AIC vs BIC)
#'    * `seed` seed for generating the data
#'    * `var.lev` number of variables and levels of the variables.
#'    * `brange` range of beta parameters
#'    * `lrange` sum of `lambdarange`
#'    * `n` sample size
#'    * `bsat` number of parameters in saturated model
#'    * `bstart` number of parameters in starting model
#'    * `edf` effective degrees of freedom
#'    * `logl` maximized log-likelihood
#'    * `dev` deviance
#'    * `icvalue` value of AIC/BIC
#'    * `pval` p-value of the fitted model
#'    * `par0` number of beta parameters with value 0
#'    * `to_0` number of parameters in `par0` estimated as 0
#'    * `budget` sum of absolute values of the estimated beta parameters
#'    * `lambdanr` number of `lambda` in the sequence of the tested lambda values
#' * `bhats` data frame with the sampled and estimated beta parameters.
#'
#'
#' @examples
#' fit(simdat())
#'
#' @details
#' The starting model is a hierarchical model excluding the k-factor interactions for which
#' one or more of the constituting (k - 1)-interactions have a non-zero sufficient statistic.
#'
#' The models with an L0 penalty include `L0a` and `Lob`. For the former the BIC is computed as
#' \eqn{-2*\ell(\beta)+edf*\log(n)}, where `edf` are the effective degrees of freedom of
#' a ridge penalty, while for the latter the BIC is computed as \eqn{-2*\ell(\beta)+df*\log(n)},
#' where `df` are the number parameters with an absolute value greater than 1e-7.
#'
#' @references Frommlet, F. and Nuel, G. (2016). An adaptive ridge regression procedure for
#' *L0* regularization. *PLOS ONE*, **11** (2), <doi:10.1371/journal.pone.0148620>.
#'
#' @importFrom stats model.matrix formula dpois rmultinom coef terms extractAIC glm poisson predict step logLik
#' @importFrom dplyr full_join join_by starts_with
#' @importFrom tidyr replace_na drop_na
#'
#' @export

fit <- function(object, lambdarange = c(1e-4, 1e-2), B = 50)
{
  # control parameters

  q      <-  0
  gamma  <-  2
  delta  <-  1e-5
  crit   <-  1e-5
  maxit  <-  250

  data   <- object$data
  f      <- object$fmax
  beta   <- object$beta

  D      <- model.matrix(f, data = data)
  nobs   <- data$Freq
  n      <- sum(nobs)
  y      <- nobs / n


  # initialize container objects
  lambdapath <- 10^seq(log10(lambdarange[1]), log10(lambdarange[2]), length.out = B)
  betas      <- matrix(nrow = ncol(D), ncol = 0,
                       dimnames = list(colnames(D), NULL))

  ################################
  # L0 estimates                 #
  ################################

  q2         <- q - 2
  w          <- c(0, rep(1, ncol(D) - 1))
  b          <- c(-log(nrow(D)), rep(0, ncol(D) - 1))
  betapath   <- edfpath <- loglpath <- gofpath <- NULL

  for (lambda in lambdapath)
  {
    tol  <- 1
    iter <- 1

    while(tol > crit)
    {
      mu  <- exp(D %*% b)[, 1]
      g   <- crossprod(D, y - mu) - 2 * lambda * w * b
      h   <- -crossprod(D, diag(mu)) %*% D - 2 * lambda * diag(w)
      bt  <- tryCatch(b - solve(h, g), error = function(e) return(b / 2))
      for (i in 2:length(w))
      {
        if (abs(bt[i]) <= delta) {
          w[i] <- delta^q2 * exp(log1p(abs(bt[i] / delta)^gamma) * q2 / gamma)
        } else {
          w[i] <- abs(bt[i])^q2 * exp(log1p(abs(delta / bt[i])^gamma) * q2 / gamma)
        }
      }
      tol  <- sum(abs(b - bt))
      iter <- iter + 1
      if (iter == maxit) break
      b    <- bt
    }

    # update values
    b        <- bt
    mu       <- exp(D %*% b)[, 1]

    # compute effective degrees of freedom and BIC
    tDW      <- crossprod(D, diag(mu))
    tmp      <- tryCatch(solve(tDW %*% D + 2 * lambda * diag(w)) %*% tDW %*% D,
                         error = function(e) return(length(b)))
    edf      <- sum(diag(tmp))
    logl     <- sum(dpois(nobs, n * mu, T))

    # update container objects
    betapath <- cbind(betapath, b)
    edfpath  <- c(edfpath, edf)
    loglpath <- c(loglpath, logl)
    gofpath  <- rbind(gofpath, c(2 * (edf - logl), edf * log(n) - 2 * logl))
  }
  # add log(n) to the intercepts in betapath

  betapath[1, ] <- betapath[1, ] + log(n)

  # best L0 model
  bestaic <- which.min(gofpath[, 1])
  bestbic <- which.min(gofpath[, 2])

  #################################
  # Step estimates                #
  #################################

  m0      <- suppressWarnings(glm(f, poisson, data = data))
  maic    <- suppressWarnings(step(m0, trace = 0))
  stepaic <- maic$coefficients
  mbic    <- suppressWarnings(step(maic, trace = 0, k = log(n)))
  stepbic <- mbic$coefficients


  #################################
  # Collect output                #
  #################################

  bhats   <- data.frame(pars = names(beta), beta = beta) %>%
    full_join(., data.frame(pars   = rownames(betapath),
                            L0aic  = round(betapath[, bestaic], 7),
                            L0bic  = round(betapath[, bestbic], 7)),
              by = join_by(pars)) %>%
    full_join(., data.frame(pars = names(stepaic), stepaic = round(stepaic, 7)),
              by = join_by(pars)) %>%
    full_join(., data.frame(pars = names(stepbic), stepbic = round(stepbic, 7)),
              by = join_by(pars)) %>%
    drop_na(starts_with("L")) %>%
    replace_na(list(stepaic = 0, stepbic = 0))

  dm  <- data.matrix(bhats[, 3:6])
  m   <- exp(D %*% dm)
  dev <- 2 * colSums(nobs * log(nobs / m), na.rm = T)

  to_0     <-  c(sum(bhats[, "beta"] == 0 & bhats[, "L0aic"] == 0),
                 sum(bhats[, "beta"] == 0 & bhats[, "L0bic"] == 0),
                 sum(bhats[, "beta"] == 0 & bhats[, "stepaic"] == 0),
                 sum(bhats[, "beta"] == 0 & bhats[, "stepaic"] == 0))

  edf       <- c(edfpath[c(bestaic, bestbic)], length(beta) - c(maic$df.res, mbic$df.res))

  # get results for lambda with lowest BIC
  out <- data.frame(method   = c("L0", "L0", "step", "step"),
                    ic       = c("aic", "bic", "aic", "bic"),
                    seed     = object$seed,
                    var.lev  = length(all.vars(f)) - 1 + length(levels(data$A)) / 10,
                    brange   = max(beta[-1]) - min(beta[-1]),
                    lrange   = sum(lambdarange),
                    n        = n,
                    bsat     = length(beta),
                    bstart   = nrow(bhats),
                    edf      = round(edf, 1),
                    logl     = round(c(loglpath[c(bestaic, bestbic)], logLik(maic)[1], logLik(mbic)[1]), 1),
                    dev      = round(dev, 1),
                    icvalue  = round(c(gofpath[bestaic, 1], gofpath[bestbic, 2], maic$aic, extractAIC(mbic, k = log(n))[2]), 1),
                    pval     = round(pchisq(dev, nrow(bhats) - edf, lower.tail = F), 3),
                    par0     = sum(bhats$beta == 0),
                    to_0     = to_0,
                    budget   = round(colSums(abs(bhats[3:6])), 1),
                    nrlambda = c(bestaic, bestbic, NA, NA)
  )
  list(out = out, bhats = bhats)
}



