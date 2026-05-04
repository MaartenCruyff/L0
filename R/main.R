#' Fit log linear models with step and L0 penalty
#'
#' @description
#' Fit log linear models with step and L0 penalty
#'
#' @param object object made with \code{\link{simdat}}.
#' @param lambdarange range of the lambda values for the L0 penalties.
#' @param B number of lambda values in `lambdarange`.
#'
#' @return
#' A list with elements
#' * `out` data frame with per method the number of non-zero parameter estimates `df`,
#' the number of effective degrees of freedom `def`, the BIC `bic`, the number of
#' zero beta parameters `par0` in the starting model, `est0` the number of which have
#' a zero estimate, `lambda` the value yielding the best BIC, and `nrlambda` its number in
#' the lambda sequence.
#' * `bhats` data frame with the sampled and estimated beta parameters.
#'
#'
#' @examples
#' fit(simdat())
#'
#' @details
#' The starting model is a hierachical model excluding the k-factor interactions for which
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
#' @importFrom stats model.matrix formula dpois rmultinom coef terms extractAIC
#' glm poisson predict step logLik
#' @importFrom dplyr full_join join_by
#' @importFrom tidyr replace_na drop_na
#'
#' @export

fit <- function(object, lambdarange = c(1e-5, 1e-2), B = 50)
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
  betapath   <- edfpath <- bicpath <- loglpath <- NULL

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
    df0      <- sum(round(b, 7) != 0)
    logl     <- sum(dpois(nobs, n * mu, T))

    # update container objects
    betapath <- cbind(betapath, b)
    edfpath  <- rbind(edfpath, c(edf, df0))
    loglpath <- c(loglpath, logl)
    bicpath  <- rbind(bicpath, c(edf * log(n) - 2 * logl, df0 *log(n) - 2 * logl))
  }

  # best L0 model
  bestedf <- which.min(bicpath[, 1])
  bestdf0 <- which.min(bicpath[, 2])

  #################################
  # Step estimates                #
  #################################

  mg      <- glm(f, poisson, data = data)
  mgstep  <- step(mg, trace = 0, k = log(n))
  mgbic   <- extractAIC(mgstep, k = log(n))
  stepfit <- mgstep$coefficients


  #################################
  # Collect output                #
  #################################

  bhats   <- data.frame(pars = names(beta), beta = beta) %>%
    full_join(., data.frame(pars = rownames(betapath),
                            L0a  = round(betapath[, bestedf], 7),
                            L0b  = round(betapath[, bestdf0], 7)),
              by = join_by(pars)) %>%
    full_join(., data.frame(pars = names(stepfit), step = round(stepfit, 7)),
              by = join_by(pars)) %>%
    drop_na(L0a) %>%
    replace_na(list(step = 0))

  to_0     <- c(sum(bhats[, "beta"] == 0 & bhats[, "L0a"] == 0),
                sum(bhats[, "beta"] == 0 & bhats[, "L0b"] == 0),
                sum(bhats[, "beta"] == 0 & bhats[, "step"] == 0))

  # get results for lambda with lowest BIC
  out <- data.frame(method   = c("L0a", "L0b", "step"),
                    simpar   = length(beta),
                    in_model = nrow(bhats),
                    edf      = round(c(edfpath[bestedf, 1], edfpath[bestdf0, 2], mgbic[1]), 1),
                    logl     = round(c(loglpath[bestedf], loglpath[bestdf0], logLik(mgstep)[1]), 1),
                    bic      = round(c(bicpath[bestedf, 1], bicpath[bestdf0, 2], mgbic[2]), 1),
                    par0     = sum(bhats$beta == 0),
                    est0     = to_0,
                    lambda   = c(lambdapath[bestedf], lambdapath[bestdf0], NA),
                    nrlambda = c(bestedf, bestdf0, NA)
  )
  list(out = out, bhats = bhats)
}




