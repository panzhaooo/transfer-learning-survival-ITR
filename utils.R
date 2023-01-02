##################
## data generation
simData <- function(N = 50000, m = 2000,
                    random_assignment = FALSE) {
  cov.mat <- matrix(c(1, 0, 0.2,
                      0, 1, 0,
                      0.2, 0, 1), 3, 3)
  # target sample
  Xtarget_population <- mvrnorm(N, mu = rep(0, 3), Sigma = cov.mat)
  Xtarget <- Xtarget_population[sample(N, m), ]
  rm(Xtarget_population)
  Xtarget[Xtarget > 4] <- 4
  Xtarget[Xtarget < -4] <- -4
  colnames(Xtarget) <- paste0("X", 1:3)
  
  # source sample
  Xtarget_population <- mvrnorm(N, mu = rep(0, 3), Sigma = cov.mat)
  Xtarget_population[Xtarget_population > 4] <- 4
  Xtarget_population[Xtarget_population < -4] <- -4
  # sampling score
  ss <- plogis(-4.5 - 0.5 * Xtarget_population[, 1] 
               - 0.5 * Xtarget_population[, 2] - 0.4 * Xtarget_population[, 3])
  Xsource <- Xtarget_population[runif(N) < ss, ]
  rm(Xtarget_population, ss)
  colnames(Xsource) <- paste0("X", 1:3)
  n_source <- nrow(Xsource)
  
  if (random_assignment) {
    ps <- 0.5 # propensity score
  } else {
    ps <- plogis(0.5 + 0.8 * Xsource[, 1] 
                 - 0.5 * Xsource[, 2] + 0 * Xsource[, 3])
  }
  A <- as.numeric(runif(n_source) < ps)
  
  u <- runif(n_source)
  lambdaT <- 1
  T0 <- log(1 - log(u) / lambdaT / exp(-2.5 - 1.5 * Xsource[, 1] 
                                       - 1 * Xsource[, 2] - 0.7 * Xsource[, 3]))
  T1 <- log(1 - log(u) / lambdaT / exp(-1 - 1 * Xsource[, 1] 
                                       - 0.9 * Xsource[, 2] - 1 * Xsource[, 3] 
                                       - 2 * Xsource[, 2]^2 
                                       + Xsource[, 1] * Xsource[, 3]))
  Ta <- T1 * A + T0 * (1 - A)
  
  u <- runif(n_source)
  lambdaC <- 0.04
  C0 <- log(1 - log(u) / lambdaC / exp(-1.6 + 0.8 * Xsource[, 1] 
                                       - 1.1 * Xsource[, 2] - 0.7 * Xsource[, 3]))
  C1 <- log(1 - log(u) / lambdaC / exp(-1.8 - 0.8 * Xsource[, 1] 
                                       - 1.7 * Xsource[, 2] - 1.4 * Xsource[, 3]))
  Ca <- C1 * A + C0 * (1 - A)
  
  U <- pmin(Ta, Ca)
  status <- as.numeric(Ta <= Ca)
  
  out <- list(Xtarget = Xtarget, Xsource = Xsource, 
              A = A, U = U, status = status)
  return(out)
}


##############################
## propensity score estimation
compute_ps <- function(X, A, 
                       model = "LR", misspecification = FALSE,
                       crossfit = NULL) {
  # model: logistic regression (LR) or random forest (RF)
  if (model == "LR") {
    if (misspecification) {
      glmA <- glm(A ~ exp(X[, 3]), family = binomial(link = "logit"))
    } else {
      glmA <- glm(A ~ X[, 1:2], family = binomial(link = "logit"))
    }
    out <- glmA$fitted.values
  } else if (model == "RF") {
    A <- as.factor(A)
    if (is.null(crossfit)) {
      grfA <- probability_forest(X, A, honesty = FALSE, num.trees = 2000,
                                 sample.fraction = 0.8, ci.group.size = 1)
      out <- c(grfA$predictions[, 2])
    } else {
      out <- rep(NA, length(A))
      for (fold in unique(crossfit)) {
        grfA <- probability_forest(X[crossfit != fold, ], A[crossfit != fold],
                                   honesty = FALSE, num.trees = 2000,
                                   compute.oob.predictions = FALSE,
                                   sample.fraction = 0.9, ci.group.size = 1)
        out[crossfit == fold] <- c(predict(grfA, X[crossfit == fold, ])$predictions[, 2])
      }
    }
  } else {
    stop("model must be set to LR or RF")
  }
  return(out)
}


################################
## survival functions estimation
survivalFit <- function(U, A, X, status, newX,
                        model = "Cox", 
                        event_misspecification = FALSE,
                        censor_misspecification = FALSE,
                        crossfit = NULL) {
  e <- length(U)
  ld <- lower.tri(matrix(1, e, e), diag = TRUE) * 1
  ud <- matrix(1, e, e) - ld
  # model: Cox proportional hazards (Cox) or random forest (RF)
  if (model == "Cox") {
    # code from Bai et al. (2017) doi: 10.1007/s10985-016-9376-x
    e <- length(U)
    # compute Cox model coefficients
    U0 <- U[A == 0]
    status0 <- status[A == 0]
    Xa0 <- X[A == 0, ]
    U1 <- U[A == 1]
    status1 <- status[A == 1]
    Xa1 <- X[A == 1, ]
    if (event_misspecification) {
      coefT0 <- coxph(Surv(U0, status0) ~ exp(Xa0[, 1:3]))$coef
      coefT1 <- coxph(Surv(U1, status1) ~ exp(Xa1[, 1:3]))$coef
    } else {
      coefT0 <- coxph(Surv(U0, status0) ~ Xa0[, 1:3])$coef
      coefT1 <- coxph(Surv(U1, status1) ~ Xa1[, 1] + Xa1[, 2] + Xa1[, 3] 
                      + I(Xa1[, 2]^2) + Xa1[, 1]:Xa1[, 3])$coef
    }
    coefT0[is.na(coefT0)] <- 0
    coefT1[is.na(coefT1)] <- 0
    if (censor_misspecification) {
      coefC0 <- coxph(Surv(U0, 1 - status0) ~ exp(Xa0[, 1:3]))$coef
      coefC1 <- coxph(Surv(U1, 1 - status1) ~ exp(Xa1[, 1:3]))$coef
    } else {
      coefC0 <- coxph(Surv(U0, 1 - status0) ~ Xa0[, 1:3])$coef
      coefC1 <- coxph(Surv(U1, 1 - status1) ~ Xa1[, 1:3])$coef
    }
    coefC0[is.na(coefC0)] <- 0
    coefC1[is.na(coefC1)] <- 0
    
    if (event_misspecification) {
      expT0 <- exp(exp(X[, 1:3]) %*% coefT0)
      expT1 <- exp(exp(X[, 1:3]) %*% coefT1)
    } else {
      expT0 <- exp(X[, 1:3] %*% coefT0)
      expT1 <- exp(cbind(X[, 1:3], X[, 2]^2, X[, 1] * X[, 3]) %*% coefT1)
    }
    if (censor_misspecification) {
      expC0 <- exp(exp(X[, 1:3]) %*% coefC0)
      expC1 <- exp(exp(X[, 1:3]) %*% coefC1)
    } else {
      expC0 <- exp(X[, 1:3] %*% coefC0)
      expC1 <- exp(X[, 1:3] %*% coefC1)
    }
    
    # compute proportional hazards
    sT0 <- cumsum((1 - A[e:1]) * expT0[e:1])
    sT1 <- cumsum(A[e:1] * expT1[e:1])
    dls0 <- (1 - A) * status / sT0[e:1]
    dls0[is.nan(dls0)] <- 0
    dls1 <- A * status / sT1[e:1]
    dls1[is.nan(dls1)] <- 0
    
    sC0 <- cumsum((1 - A[e:1]) * expC0[e:1])
    sC1 <- cumsum(A[e:1] * expC1[e:1])
    dlc0 <- (1 - A) * (1 - status) / sC0[e:1]
    dlc0[is.nan(dlc0)] <- 0
    dlc1 <- A * (1 - status) / sC1[e:1]
    dlc1[is.nan(dlc1)] <- 0
    
    # compute survival functions
    sss0 <- matrix(0, e, e)
    ssc0 <- matrix(0, e, e)
    sss1 <- matrix(0, e, e)
    ssc1 <- matrix(0, e, e)
    
    sss0[, 1] <- dls0[1] * expT0
    ssc0[, 1] <- dlc0[1] * expC0
    
    ots0 <- outer(c(expT0), dls0, "*")
    ots0[, 1] <- rep(0, e)
    sss0 <- outer(sss0[, 1], rep(1, e)) + t(apply(ots0, 1, cumsum))
    
    otc0 <- outer(c(expC0), dlc0, "*")
    otc0[, 1] <- rep(0, e)
    ssc0 <- outer(ssc0[, 1], rep(1, e)) + t(apply(otc0, 1, cumsum))
    
    sss1[, 1] <- dls1[1] * expT1
    ssc1[, 1] <- dlc1[1] * expC1
    
    ots1 <- outer(c(expT1), dls1, "*")
    ots1[, 1] <- rep(0, e)
    sss1 <- outer(sss1[, 1], rep(1, e)) + t(apply(ots1, 1, cumsum))
    
    otc1 <- outer(c(expC1), dlc1, "*")
    otc1[, 1] <- rep(0, e)
    ssc1 <- outer(ssc1[, 1], rep(1, e)) + t(apply(otc1, 1, cumsum))
    
    S0 <- exp(-sss0)
    Sc0 <- exp(-ssc0)
    S1 <- exp(-sss1)
    Sc1 <- exp(-ssc1)
    
    # compute survival functions for newX
    if (event_misspecification) {
      expT0_new <- exp(exp(newX[, 1:3]) %*% coefT0)
      expT1_new <- exp(exp(newX[, 1:3]) %*% coefT1)
    } else {
      expT0_new <- exp(newX[, 1:3] %*% coefT0)
      expT1_new <- exp(cbind(newX[, 1:3], newX[, 2]^2, newX[, 1] * newX[, 3]) %*% coefT1)
    }
    e_new <- nrow(newX)
    sss0_new <- matrix(0, e_new, e)
    sss1_new <- matrix(0, e_new, e)
    
    sss0_new[, 1] <- dls0[1] * expT0_new
    ots0_new <- outer(c(expT0_new), dls0, "*")
    ots0_new[, 1] <- rep(0, e_new)
    sss0_new <- outer(sss0_new[, 1], rep(1, e)) + t(apply(ots0_new, 1, cumsum))
    
    sss1_new[, 1] <- dls1[1] * expT1_new
    ots1_new <- outer(c(expT1_new), dls1, "*")
    ots1_new[, 1] <- rep(0, e_new)
    sss1_new <- outer(sss1_new[, 1], rep(1, e)) + t(apply(ots1_new, 1, cumsum))
    
    S0_new <- exp(-sss0_new)
    S1_new <- exp(-sss1_new)
    
    # compute martingale integral
    q0 <- matrix(0, e, e)
    otr0 <- outer(c(expC0), dlc0) / (S0 * Sc0)
    otr0[is.na(otr0)] <- 0
    q0 <- t(apply(otr0, 1, cumsum))
    q0 <- q0 * ld + diag(q0) * ud
    q0[is.na(q0)] <- 0
    
    q1 <- matrix(0, e, e)
    otr1 <- outer(c(expC1), dlc1) / (S1 * Sc1)
    otr1[is.na(otr1)] <- 0
    q1 <- t(apply(otr1, 1, cumsum))
    q1 <- q1 * ld + diag(q1) * ud
    q1[is.na(q1)] <- 0
  } else if (model == "RF") {
    forestT <- survival_forest(X = data.frame(A = A, X), Y = U, D = status,
                               failure.times = U, num.trees = 2000, 
                               sample.fraction = 0.9, prediction.type = "Nelson-Aalen",
                               honesty = FALSE, compute.oob.predictions = FALSE)
    S1_new <- predict(forestT, newdata = data.frame(A = 1, newX), failure.times = U)$predictions
    S0_new <- predict(forestT, newdata = data.frame(A = 0, newX), failure.times = U)$predictions
    if (is.null(crossfit)) {
      S1 <- predict(forestT, newdata = data.frame(A = 1, X), failure.times = U)$predictions
      S0 <- predict(forestT, newdata = data.frame(A = 0, X), failure.times = U)$predictions
      forestC <- survival_forest(X = data.frame(A = A, X), Y = U, D = 1 - status,
                                 failure.times = U, num.trees = 2000, 
                                 sample.fraction = 0.9, prediction.type = "Nelson-Aalen",
                                 honesty = FALSE, compute.oob.predictions = FALSE)
      Sc1 <- predict(forestC, newdata = data.frame(A = 1, X), failure.times = U)$predictions
      Sc0 <- predict(forestC, newdata = data.frame(A = 0, X), failure.times = U)$predictions
    } else {
      S1 <- matrix(0, e, e)
      S0 <- matrix(0, e, e)
      Sc1 <-matrix(0, e, e)
      Sc0 <- matrix(0, e, e)
      for (fold in unique(crossfit)) {
        forestT <- survival_forest(X = data.frame(A = A[crossfit != fold],
                                                  X[crossfit != fold, ]), 
                                   Y = U[crossfit != fold], 
                                   D = status[crossfit != fold],
                                   failure.times = U, num.trees = 2000, 
                                   sample.fraction = 0.8, prediction.type = "Nelson-Aalen",
                                   honesty = FALSE, compute.oob.predictions = FALSE)
        S1[crossfit == fold, ] <- predict(forestT, failure.times = U,
                                          newdata = data.frame(A = 1, X[crossfit == fold, ]))$predictions
        S0[crossfit == fold, ] <- predict(forestT, failure.times = U,
                                          newdata = data.frame(A = 0, X[crossfit == fold, ]))$predictions
        forestC <- survival_forest(X = data.frame(A = A[crossfit != fold],
                                                  X[crossfit != fold, ]), 
                                   Y = U[crossfit != fold], 
                                   D = 1 - status[crossfit != fold],
                                   failure.times = U, num.trees = 2000, 
                                   sample.fraction = 0.8, prediction.type = "Nelson-Aalen",
                                   honesty = FALSE, compute.oob.predictions = FALSE)
        Sc1[crossfit == fold, ] <- predict(forestC, failure.times = U,
                                           newdata = data.frame(A = 1, X[crossfit == fold, ]))$predictions
        Sc0[crossfit == fold, ] <- predict(forestC, failure.times = U,
                                           newdata = data.frame(A = 0, X[crossfit == fold, ]))$predictions
      }
    }
    dlc0 <- t(apply(cbind(0, -log(Sc0)), 1, diff))
    dlc1 <- t(apply(cbind(0, -log(Sc1)), 1, diff))
    
    otr0 <- dlc0 / (S0 * Sc0)
    otr0[is.na(otr0)] <- 0
    q0 <- t(apply(otr0, 1, cumsum))
    q0 <- q0 * ld + diag(q0) * ud
    q0[is.na(q0)] <- 0
    
    otr1 <- dlc1 / (S1 * Sc1)
    otr1[is.na(otr1)] <- 0
    q1 <- t(apply(otr1, 1, cumsum))
    q1 <- q1 * ld + diag(q1) * ud
    q1[is.na(q1)] <- 0
  } else {
    stop("model must be set to Cox or RF")
  }
  return(list(S0 = S0, S1 = S1,
              Sc0 = Sc0, Sc1 = Sc1,
              S0_new = S0_new, S1_new = S1_new,
              q0 = q0, q1 = q1))
}


######################
## calibration weights
lamFun <- function(lam, moments, moments.bar) {
  qi <- (exp(moments %*% lam) / sum(exp(moments %*% lam)))[, 1]
  colSums(qi * moments) - moments.bar
}

compute_calibration_weights <- function(Xtarget, Xsource,
                                        model = "d1",
                                        misspecification = FALSE,
                                        crossfit = NULL) {
  if (model == "d1") {
    if (misspecification) {
      Xtarget <- exp(Xtarget[, 1, drop = FALSE])
      Xsource <- exp(Xsource[, 1, drop = FALSE])
    }
    Mbar <- colMeans(Xtarget)
    lam.hat <- searchZeros(matrix(rnorm(length(Mbar) * 5, 0, 0.5), nrow = 5), 
                           lamFun, moments = Xsource, moments.bar = Mbar)$x[1, ]
    out <- c(exp(Xsource %*% lam.hat) / sum(exp(Xsource %*% lam.hat)))
  } else if (model == "d2") {
    Mtarget <- cbind(Xtarget, Xtarget^2)
    Msource <- cbind(Xsource, Xsource^2)
    if (is.null(crossfit)) {
      Mbar <- colMeans(Mtarget)
      lam.hat <- searchZeros(matrix(rnorm(length(Mbar) * 5, 0, 0.5), nrow = 5), 
                             lamFun, moments = Msource, moments.bar = Mbar)$x[1, ]
      out <- c(exp(Msource %*% lam.hat) / sum(exp(Msource %*% lam.hat)))
    } else {
      out <- rep(NA, nrow(Xsource))
      V <- length(unique(crossfit))
      cf_target <- sample(unique(crossfit), nrow(Xtarget), replace = TRUE)
      for (fold in unique(crossfit)) {
        Mbar <- colMeans(Mtarget[cf_target != fold, ])
        lam.hat <- searchZeros(matrix(rnorm(length(Mbar) * 5, 0, 0.5), nrow = 5), 
                               lamFun, moments = Msource[crossfit != fold, ], moments.bar = Mbar)$x[1, ]
        out[crossfit == fold] <- c(exp(Msource[crossfit == fold, ] %*% lam.hat) / sum(exp(Msource[crossfit == fold, ] %*% lam.hat)) 
                                   * sum(crossfit == fold) / length(crossfit))
      }
    }
  } else {
    stop("model must be set to d1 or d2")
  }
  return(out)
}


###############
## compute RMST
compute_RMST <- function(S, U, L) {
  Lidx <- findInterval(L, U)
  time_diff <- diff(c(0, U))
  RMST <- apply(cbind(1, S), 1, function(x) sum(time_diff[1:Lidx] * x[1:Lidx]) 
                + (L - U[Lidx]) * x[Lidx + 1])
  return(RMST)
}


################################
## optimal linear ITR estimation
value_estimators <- function(Xsource, A, U, status, Xtarget,
                             event_misspecification = FALSE,
                             sampling_misspecification = FALSE,
                             propensity_misspecification = FALSE, 
                             censor_misspecification = FALSE,
                             model_propensity = "LR",
                             model_survival = "Cox",
                             model_calibration = "d1",
                             crossfitV = NULL) {
  Uorder <- order(U)
  U <- U[Uorder]
  Xsource <- Xsource[Uorder, ]
  status <- status[Uorder]
  A <- A[Uorder]
  
  if (is.null(crossfitV)) {
    crossfit <- NULL
  } else {
    crossfit <- sample(1:crossfitV, length(A), replace = TRUE)
  }
  
  # nuisance parameters estimation
  ps <- compute_ps(Xsource, A, model_propensity, 
                   propensity_misspecification, crossfit)
  psA <- A / ps + (1 - A) / (1 - ps)
  
  sF <- survivalFit(U, A, Xsource, status, Xtarget,
                    model_survival,
                    event_misspecification,
                    censor_misspecification,
                    crossfit)
  Sa <- A * sF$S1 + (1 - A) * sF$S0
  Sca <- A * sF$Sc1 + (1 - A) * sF$Sc0
  
  q <- compute_calibration_weights(Xtarget, Xsource, 
                                   model_calibration, 
                                   sampling_misspecification,
                                   crossfit)
  
  n <- nrow(Xsource)
  m <- nrow(Xtarget)
  ld <- lower.tri(matrix(1, n, n), diag = TRUE) * 1
  ud <- matrix(1, n, n) - ld
  
  # RMST, max time L
  UL <- pmin(U, L)
  statusL <- as.numeric(U > L) + status * as.numeric(U <= L)
  sca <- Sca[cbind(1:n, findInterval(UL, U))] # at time Ui
  
  # naive
  RMSTnaive <- psA * statusL / sca * UL
  Vhat <- function(d) {
    dA <- as.numeric(cbind(1, Xsource) %*% d > 0)
    mean(as.numeric(A == dA) * RMSTnaive)
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vnaive <- Vopt$value
  d_naive <- Vopt$par
  d_naive <- d_naive / sqrt(sum(d_naive^2))
  
  # IPSW
  V <- c(rep(1, n), rep(0, m)) # indicator of source sample
  # estimate P(V = 1 | X)
  dataV <- rbind(Xsource, Xtarget) 
  dataV <- as.data.frame(cbind(V, dataV))
  if (sampling_misspecification) {
    glmV <- glm(V ~ X1, data = dataV, family = binomial(link = "logit"))
  } else {
    glmV <- glm(V ~ X1 + X2 + X3, data = dataV, family = binomial(link = "logit"))
  }
  p <- glmV$fitted.values[1:n]
  odds <- (1 - p) / p
  
  RMSTipsw <- odds * psA * statusL / sca * UL
  Vhat <- function(d) {
    dA <- as.numeric(cbind(1, Xsource) %*% d > 0)
    sum(as.numeric(A == dA) * RMSTipsw) / m
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vipsw <- Vopt$value
  d_ipsw <- Vopt$par
  d_ipsw <- d_ipsw / sqrt(sum(d_ipsw^2))
  
  # CW_IPW
  RMSTcwipw <- n * q * psA * statusL / sca * UL
  Vhat <- function(d) {
    dA <- as.numeric(cbind(1, Xsource) %*% d > 0)
    mean(as.numeric(A == dA) * RMSTcwipw)
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vcwipw <- Vopt$value
  d_cwipw <- Vopt$par
  d_cwipw <- d_cwipw / sqrt(sum(d_cwipw^2))
  
  # CW_OR
  RMSTcwor0 <- q * compute_RMST(sF$S0, U, L)
  RMSTcwor1 <- q * compute_RMST(sF$S1, U, L)
  Vhat <- function(d) {
    dA <- as.numeric(cbind(1, Xsource) %*% d > 0)
    sum(dA * RMSTcwor1 + (1 - dA) * RMSTcwor0)
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vcwor <- Vopt$value
  d_cwor <- Vopt$par
  d_cwor <- d_cwor / sqrt(sum(d_cwor^2))
  
  # OR target
  RMSTort0 <- compute_RMST(sF$S0_new, U, L)
  RMSTort1 <- compute_RMST(sF$S1_new, U, L)
  Vhat <- function(d) {
    dA <- as.numeric(cbind(1, Xtarget) %*% d > 0)
    mean(dA * RMSTort1 + (1 - dA) * RMSTort0)
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vort <- Vopt$value
  d_ort <- Vopt$par
  d_ort <- d_ort / sqrt(sum(d_ort^2))
  
  # ACW
  Saug <- psA * (ld / Sca - Sa + (1 - status) * Sa / (diag(Sa) * diag(Sca)) * ud
                 - Sa * (A * sF$q1 + (1 - A) * sF$q0)) 
  Saug[is.na(Saug)] <- 0
  Saug[is.infinite(Saug)] <- 0
  Saug[is.nan(Saug)] <- 0
  RMSTaug <- compute_RMST(Saug, U, L)
  Vhat <- function(d) {
    dA1 <- as.numeric(cbind(1, Xtarget) %*% d > 0)
    dA2 <- as.numeric(cbind(1, Xsource) %*% d > 0)
    mean(dA1 * RMSTort1 + (1 - dA1) * RMSTort0) + sum(as.numeric(A == dA2) * q * RMSTaug)
  }
  Vopt <- genoud(Vhat, default.domains = 3000, nvars = 4,
                 pop.size = 6000, max = TRUE, print.level = 0,
                 wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
                 P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
                 P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
                 starting.values = NULL, hard.generation.limit = FALSE,
                 solution.tolerance = 1e-04, optim.method = "Nelder-Mead")
  Vacw <- Vopt$value
  d_acw <- Vopt$par
  d_acw <- d_acw / sqrt(sum(d_acw^2))
  
  return(list(Vnaive = Vnaive, d_naive = d_naive,
              Vipsw = Vipsw, d_ipsw = d_ipsw,
              Vcwipw = Vcwipw, d_cwipw = d_cwipw,
              Vcwor = Vcwor, d_cwor = d_cwor,
              Vort = Vort, d_ort = d_ort,
              Vacw = Vacw, d_acw = d_acw))
}





