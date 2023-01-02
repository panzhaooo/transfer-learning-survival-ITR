library(MASS)
library(rgenoud)
library(survival)
library(nleqslv)
library(doParallel)
library(grf)

##########################
## (semi)parametric models

# simulation setup
N <- 2e5
m <- 8000
random_treatment <- FALSE

L <- 4

model_propensity <- "LR"
model_survival <- "Cox"
model_calibration <- "d1"
crossfitV <- NULL

n.MonteCarlo <- 350
B <- 200

# generate a target population for evaluation
set.seed(1997)
largeN <- 1e5
largeXtarget <- mvrnorm(largeN, mu = rep(0, 3), Sigma = matrix(c(1, 0, 0.2,
                                                                 0, 1, 0,
                                                                 0.2, 0, 1), 3, 3))
largeXtarget[largeXtarget > 4] <- 4
largeXtarget[largeXtarget < -4] <- -4

u <- runif(largeN)
lambdaT <- 1
largeT0 <- log(1 - log(u) / lambdaT / exp(-2.5 - 1.5 * largeXtarget[, 1] - largeXtarget[, 2] - 0.7 * largeXtarget[, 3]))
largeT1 <- log(1 - log(u) / lambdaT / exp(-1 - largeXtarget[, 1] 
                                          - 0.9 * largeXtarget[, 2] - 1 * largeXtarget[, 3] 
                                          - 2 * largeXtarget[, 2]^2 + largeXtarget[, 1] * largeXtarget[, 3]))

# approximate true optimal linear ITR
largeT0L <- pmin(largeT0, L)
largeT1L <- pmin(largeT1, L)
Vd <- function(d) {
  Ad <- as.numeric(cbind(1, largeXtarget) %*% d > 0)
  mean(Ad * largeT1L + (1 - Ad) * largeT0L)
}

# parallel computation
nworkers <- detectCores()
cl <- makePSOCKcluster(nworkers - 1)
clusterExport(cl, c("largeXtarget", "largeT1L", "largeT0L"))

Vmax <- genoud(Vd, default.domains = 3000, nvars = 4, 
               pop.size = 10000, max = TRUE, print.level = 0,
               wait.generations = 8, gradient.check = FALSE, BFGS = FALSE,
               P1 = 50, P2 = 50, P3 = 10, P4 = 50, P5 = 50,
               P6 = 50, P7 = 50, P8 = 50, P9 = 0, P9mix = NULL,
               starting.values = NULL, hard.generation.limit = FALSE,
               solution.tolerance = 1e-04, optim.method = "Nelder-Mead",
               cluster = cl)
stopCluster(cl)
cat(0, as.character(Sys.time()), nworkers, "\n")
# end
d_opt <- Vmax$par
d_opt <- d_opt / sqrt(sum(d_opt ^ 2))
Aopt <- as.numeric(cbind(1, largeXtarget) %*% d_opt > 0)
Vmax <- Vmax$value

pcd <- function(d) {
  Ad <- as.numeric(cbind(1, largeXtarget) %*% d > 0)
  mean(Ad == Aopt)
}

#########
## CASE 1

Vd_naive <- Vd_ipsw <- Vd_cwipw <- Vd_cwor <- Vd_ort <- Vd_acw <- rep(0, n.MonteCarlo)
pcd_naive <- pcd_ipsw <- pcd_cwipw <- pcd_cwor <- pcd_ort <- pcd_acw <- rep(0, n.MonteCarlo)
Vhat_naive <- Vhat_ipsw <- Vhat_cwipw <- Vhat_cwor <- Vhat_ort <- Vhat_acw <- rep(0, n.MonteCarlo)

res.naive <- matrix(0, n.MonteCarlo, 4)
colnames(res.naive) <- c("bias", "se", "coverageW", "coverageQ")
res.ipsw <- res.cw.ipw <- res.cw.or <- res.ort <- res.acw <- res.naive

event_misspecification <- FALSE
sampling_misspecification <- FALSE
propensity_misspecification <- FALSE
censor_misspecification <- FALSE

set.seed(6)
for (i in 1:n.MonteCarlo) {
  dat <- simData(N, m, random_treatment)
  estimators <- value_estimators(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
                                 event_misspecification,
                                 sampling_misspecification,
                                 propensity_misspecification, 
                                 censor_misspecification,
                                 model_propensity,
                                 model_survival,
                                 model_calibration,
                                 crossfitV)
  # true value
  Vd_naive[i] <- Vd(estimators$d_naive)
  Vd_ipsw[i] <- Vd(estimators$d_ipsw)
  Vd_cwipw[i] <- Vd(estimators$d_cwipw)
  Vd_cwor[i] <- Vd(estimators$d_cwor)
  Vd_ort[i] <- Vd(estimators$d_ort)
  Vd_acw[i] <- Vd(estimators$d_acw)
  
  # proportion of correct decision
  pcd_naive[i] <- pcd(estimators$d_naive)
  pcd_ipsw[i] <- pcd(estimators$d_ipsw)
  pcd_cwipw[i] <- pcd(estimators$d_cwipw)
  pcd_cwor[i] <- pcd(estimators$d_cwor)
  pcd_ort[i] <- pcd(estimators$d_ort)
  pcd_acw[i] <- pcd(estimators$d_acw)
  
  Vhat_naive[i] <- estimators$Vnaive
  Vhat_ipsw[i] <- estimators$Vipsw
  Vhat_cwipw[i] <- estimators$Vcwipw
  Vhat_cwor[i] <- estimators$Vcwor
  Vhat_ort[i] <- estimators$Vort
  Vhat_acw[i] <- estimators$Vacw
  
  # bootstrap variance estimation
  # parallel computation
  # n.workers <- detectCores()
  cl <- makePSOCKcluster(48)
  clusterExport(cl, c("event_misspecification", "sampling_misspecification",
                      "propensity_misspecification", "censor_misspecification",
                      "dat",
                      "model_propensity", "model_survival",
                      "model_calibration", "crossfitV",
                      "simData", "compute_ps", "survivalFit",
                      "compute_RMST", "lamFun", "compute_calibration_weights",
                      "value_estimators", "Vd", "pcd", "L",
                      "largeXtarget", "largeT1L", "largeT0L", "Aopt"))
  registerDoParallel(cl)
  
  bootV <- foreach(b = 1:B, .combine = 'cbind',
                   .packages = c("rgenoud", "survival", "nleqslv")) %dopar% {
                     n <- nrow(dat$Xsource)
                     m <- nrow(dat$Xtarget)
                     idx_b <- sample(n, replace = TRUE)
                     estimator_b <- value_estimators(dat$Xsource[idx_b, ], dat$A[idx_b], dat$U[idx_b], dat$status[idx_b], 
                                                     dat$Xtarget[sample(m, replace = TRUE), ],
                                                     event_misspecification,
                                                     sampling_misspecification,
                                                     propensity_misspecification, 
                                                     censor_misspecification,
                                                     model_propensity,
                                                     model_survival,
                                                     model_calibration,
                                                     crossfitV)
                     c(estimator_b$Vnaive, estimator_b$Vipsw, estimator_b$Vcwipw,
                       estimator_b$Vcwor, estimator_b$Vort, estimator_b$Vacw)
                   }
  stopCluster(cl)
  # end
  sd.naive <- sd(bootV[1, ])
  sd.ipsw <- sd(bootV[2, ])
  sd.cwipw <- sd(bootV[3, ])
  sd.cwor <- sd(bootV[4, ])
  sd.ort <- sd(bootV[5, ])
  sd.acw <- sd(bootV[6, ])
  
  # results summary
  res.naive[i, ] <- c(Vhat_naive[i] - Vmax, sd.naive,
                      as.numeric(Vmax >= (Vhat_naive[i] - qnorm(0.975) * sd.naive)
                                 & Vmax <= (Vhat_naive[i] + qnorm(0.975) * sd.naive)),
                      as.numeric(Vmax >= quantile(bootV[1, ], 0.025)
                                 & Vmax <= quantile(bootV[1, ], 0.975)))
  res.ipsw[i, ] <- c(Vhat_ipsw[i] - Vmax, sd.ipsw,
                     as.numeric(Vmax >= (Vhat_ipsw[i] - qnorm(0.975) * sd.ipsw)
                                & Vmax <= (Vhat_ipsw[i] + qnorm(0.975) * sd.ipsw)),
                     as.numeric(Vmax >= quantile(bootV[2, ], 0.025)
                                & Vmax <= quantile(bootV[2, ], 0.975)))
  res.cw.ipw[i, ] <- c(Vhat_cwipw[i] - Vmax, sd.cwipw,
                       as.numeric(Vmax >= (Vhat_cwipw[i] - qnorm(0.975) * sd.cwipw)
                                  & Vmax <= (Vhat_cwipw[i] + qnorm(0.975) * sd.cwipw)),
                       as.numeric(Vmax >= quantile(bootV[3, ], 0.025)
                                  & Vmax <= quantile(bootV[3, ], 0.975)))
  res.cw.or[i, ] <- c(Vhat_cwor[i] - Vmax, sd.cwor,
                      as.numeric(Vmax >= (Vhat_cwor[i] - qnorm(0.975) * sd.cwor)
                                 & Vmax <= (Vhat_cwor[i] + qnorm(0.975) * sd.cwor)),
                      as.numeric(Vmax >= quantile(bootV[4, ], 0.025)
                                 & Vmax <= quantile(bootV[4, ], 0.975)))
  res.ort[i, ] <- c(Vhat_ort[i] - Vmax, sd.ort,
                    as.numeric(Vmax >= (Vhat_ort[i] - qnorm(0.975) * sd.ort)
                               & Vmax <= (Vhat_ort[i] + qnorm(0.975) * sd.ort)),
                    as.numeric(Vmax >= quantile(bootV[5, ], 0.025)
                               & Vmax <= quantile(bootV[5, ], 0.975)))
  res.acw[i, ] <- c(Vhat_acw[i] - Vmax, sd.acw,
                    as.numeric(Vmax >= (Vhat_acw[i] - qnorm(0.975) * sd.acw)
                               & Vmax <= (Vhat_acw[i] + qnorm(0.975) * sd.acw)),
                    as.numeric(Vmax >= quantile(bootV[6, ], 0.025)
                               & Vmax <= quantile(bootV[6, ], 0.975)))
  cat(i, as.character(Sys.time()), "\n")
}

result1 <- list(Vd = list(Vd_naive = Vd_naive,
                          Vd_ipsw = Vd_ipsw,
                          Vd_cwipw = Vd_cwipw,
                          Vd_cwor = Vd_cwor,
                          Vd_ort = Vd_ort,
                          Vd_acw = Vd_acw),
                pcd = list(pcd_naive = pcd_naive,
                           pcd_ipsw = pcd_ipsw,
                           pcd_cwipw = pcd_cwipw,
                           pcd_cwor = pcd_cwor,
                           pcd_ort = pcd_ort,
                           pcd_acw = pcd_acw),
                Vhat = list(Vhat_naive = Vhat_naive,
                            Vhat_ipsw = Vhat_ipsw,
                            Vhat_cwipw = Vhat_cwipw,
                            Vhat_cwor = Vhat_cwor,
                            Vhat_ort = Vhat_ort,
                            Vhat_acw = Vhat_acw),
                res = list(res.naive = res.naive,
                           res.ipsw = res.ipsw,
                           res.cw.ipw = res.cw.ipw,
                           res.cw.or = res.cw.or,
                           res.ort = res.ort,
                           res.acw = res.acw))

save(result1, Vmax, file = "survITR1.RData")

#########
## CASE 2

Vd_naive <- Vd_ipsw <- Vd_cwipw <- Vd_cwor <- Vd_ort <- Vd_acw <- rep(0, n.MonteCarlo)
pcd_naive <- pcd_ipsw <- pcd_cwipw <- pcd_cwor <- pcd_ort <- pcd_acw <- rep(0, n.MonteCarlo)
Vhat_naive <- Vhat_ipsw <- Vhat_cwipw <- Vhat_cwor <- Vhat_ort <- Vhat_acw <- rep(0, n.MonteCarlo)

res.naive <- matrix(0, n.MonteCarlo, 4)
colnames(res.naive) <- c("bias", "se", "coverageW", "coverageQ")
res.ipsw <- res.cw.ipw <- res.cw.or <- res.ort <- res.acw <- res.naive

event_misspecification <- TRUE
sampling_misspecification <- FALSE
propensity_misspecification <- FALSE
censor_misspecification <- FALSE

set.seed(6)
for (i in 1:n.MonteCarlo) {
  dat <- simData(N, m, random_treatment)
  estimators <- value_estimators(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
                                 event_misspecification,
                                 sampling_misspecification,
                                 propensity_misspecification, 
                                 censor_misspecification,
                                 model_propensity,
                                 model_survival,
                                 model_calibration,
                                 crossfitV)
  # true value
  Vd_naive[i] <- Vd(estimators$d_naive)
  Vd_ipsw[i] <- Vd(estimators$d_ipsw)
  Vd_cwipw[i] <- Vd(estimators$d_cwipw)
  Vd_cwor[i] <- Vd(estimators$d_cwor)
  Vd_ort[i] <- Vd(estimators$d_ort)
  Vd_acw[i] <- Vd(estimators$d_acw)
  
  # proportion of correct decision
  pcd_naive[i] <- pcd(estimators$d_naive)
  pcd_ipsw[i] <- pcd(estimators$d_ipsw)
  pcd_cwipw[i] <- pcd(estimators$d_cwipw)
  pcd_cwor[i] <- pcd(estimators$d_cwor)
  pcd_ort[i] <- pcd(estimators$d_ort)
  pcd_acw[i] <- pcd(estimators$d_acw)
  
  Vhat_naive[i] <- estimators$Vnaive
  Vhat_ipsw[i] <- estimators$Vipsw
  Vhat_cwipw[i] <- estimators$Vcwipw
  Vhat_cwor[i] <- estimators$Vcwor
  Vhat_ort[i] <- estimators$Vort
  Vhat_acw[i] <- estimators$Vacw
  
  # bootstrap variance estimation
  # parallel computation
  #n.workers <- detectCores()
  cl <- makePSOCKcluster(40)
  clusterExport(cl, c("event_misspecification", "sampling_misspecification",
                      "propensity_misspecification", "censor_misspecification",
                      "dat",
                      "model_propensity", "model_survival",
                      "model_calibration", "crossfitV",
                      "simData", "compute_ps", "survivalFit",
                      "compute_RMST", "lamFun", "compute_calibration_weights",
                      "value_estimators", "Vd", "pcd", "L",
                      "largeXtarget", "largeT1L", "largeT0L", "Aopt"))
  registerDoParallel(cl)
  
  bootV <- foreach(b = 1:B, .combine = 'cbind',
                   .packages = c("rgenoud", "survival", "nleqslv")) %dopar% {
                     n <- nrow(dat$Xsource)
                     m <- nrow(dat$Xtarget)
                     idx_b <- sample(n, replace = TRUE)
                     estimator_b <- value_estimators(dat$Xsource[idx_b, ], dat$A[idx_b], dat$U[idx_b], dat$status[idx_b], 
                                                     dat$Xtarget[sample(m, replace = TRUE), ],
                                                     event_misspecification,
                                                     sampling_misspecification,
                                                     propensity_misspecification, 
                                                     censor_misspecification,
                                                     model_propensity,
                                                     model_survival,
                                                     model_calibration,
                                                     crossfitV)
                     c(estimator_b$Vnaive, estimator_b$Vipsw, estimator_b$Vcwipw,
                       estimator_b$Vcwor, estimator_b$Vort, estimator_b$Vacw)
                   }
  stopCluster(cl)
  # end
  sd.naive <- sd(bootV[1, ])
  sd.ipsw <- sd(bootV[2, ])
  sd.cwipw <- sd(bootV[3, ])
  sd.cwor <- sd(bootV[4, ])
  sd.ort <- sd(bootV[5, ])
  sd.acw <- sd(bootV[6, ])
  
  # results summary
  res.naive[i, ] <- c(Vhat_naive[i] - Vmax, sd.naive,
                      as.numeric(Vmax >= (Vhat_naive[i] - qnorm(0.975) * sd.naive)
                                 & Vmax <= (Vhat_naive[i] + qnorm(0.975) * sd.naive)),
                      as.numeric(Vmax >= quantile(bootV[1, ], 0.025)
                                 & Vmax <= quantile(bootV[1, ], 0.975)))
  res.ipsw[i, ] <- c(Vhat_ipsw[i] - Vmax, sd.ipsw,
                     as.numeric(Vmax >= (Vhat_ipsw[i] - qnorm(0.975) * sd.ipsw)
                                & Vmax <= (Vhat_ipsw[i] + qnorm(0.975) * sd.ipsw)),
                     as.numeric(Vmax >= quantile(bootV[2, ], 0.025)
                                & Vmax <= quantile(bootV[2, ], 0.975)))
  res.cw.ipw[i, ] <- c(Vhat_cwipw[i] - Vmax, sd.cwipw,
                       as.numeric(Vmax >= (Vhat_cwipw[i] - qnorm(0.975) * sd.cwipw)
                                  & Vmax <= (Vhat_cwipw[i] + qnorm(0.975) * sd.cwipw)),
                       as.numeric(Vmax >= quantile(bootV[3, ], 0.025)
                                  & Vmax <= quantile(bootV[3, ], 0.975)))
  res.cw.or[i, ] <- c(Vhat_cwor[i] - Vmax, sd.cwor,
                      as.numeric(Vmax >= (Vhat_cwor[i] - qnorm(0.975) * sd.cwor)
                                 & Vmax <= (Vhat_cwor[i] + qnorm(0.975) * sd.cwor)),
                      as.numeric(Vmax >= quantile(bootV[4, ], 0.025)
                                 & Vmax <= quantile(bootV[4, ], 0.975)))
  res.ort[i, ] <- c(Vhat_ort[i] - Vmax, sd.ort,
                    as.numeric(Vmax >= (Vhat_ort[i] - qnorm(0.975) * sd.ort)
                               & Vmax <= (Vhat_ort[i] + qnorm(0.975) * sd.ort)),
                    as.numeric(Vmax >= quantile(bootV[5, ], 0.025)
                               & Vmax <= quantile(bootV[5, ], 0.975)))
  res.acw[i, ] <- c(Vhat_acw[i] - Vmax, sd.acw,
                    as.numeric(Vmax >= (Vhat_acw[i] - qnorm(0.975) * sd.acw)
                               & Vmax <= (Vhat_acw[i] + qnorm(0.975) * sd.acw)),
                    as.numeric(Vmax >= quantile(bootV[6, ], 0.025)
                               & Vmax <= quantile(bootV[6, ], 0.975)))
  cat(i, as.character(Sys.time()), "\n")
}

result2 <- list(Vd = list(Vd_naive = Vd_naive,
                          Vd_ipsw = Vd_ipsw,
                          Vd_cwipw = Vd_cwipw,
                          Vd_cwor = Vd_cwor,
                          Vd_ort = Vd_ort,
                          Vd_acw = Vd_acw),
                pcd = list(pcd_naive = pcd_naive,
                           pcd_ipsw = pcd_ipsw,
                           pcd_cwipw = pcd_cwipw,
                           pcd_cwor = pcd_cwor,
                           pcd_ort = pcd_ort,
                           pcd_acw = pcd_acw),
                Vhat = list(Vhat_naive = Vhat_naive,
                            Vhat_ipsw = Vhat_ipsw,
                            Vhat_cwipw = Vhat_cwipw,
                            Vhat_cwor = Vhat_cwor,
                            Vhat_ort = Vhat_ort,
                            Vhat_acw = Vhat_acw),
                res = list(res.naive = res.naive,
                           res.ipsw = res.ipsw,
                           res.cw.ipw = res.cw.ipw,
                           res.cw.or = res.cw.or,
                           res.ort = res.ort,
                           res.acw = res.acw))

save(result2, file = "survITR2.RData")


#########
## CASE 3

Vd_naive <- Vd_ipsw <- Vd_cwipw <- Vd_cwor <- Vd_ort <- Vd_acw <- rep(0, n.MonteCarlo)
pcd_naive <- pcd_ipsw <- pcd_cwipw <- pcd_cwor <- pcd_ort <- pcd_acw <- rep(0, n.MonteCarlo)
Vhat_naive <- Vhat_ipsw <- Vhat_cwipw <- Vhat_cwor <- Vhat_ort <- Vhat_acw <- rep(0, n.MonteCarlo)

res.naive <- matrix(0, n.MonteCarlo, 4)
colnames(res.naive) <- c("bias", "se", "coverageW", "coverageQ")
res.ipsw <- res.cw.ipw <- res.cw.or <- res.ort <- res.acw <- res.naive

event_misspecification <- FALSE
sampling_misspecification <- TRUE
propensity_misspecification <- TRUE
censor_misspecification <- TRUE

set.seed(6)
for (i in 1:n.MonteCarlo) {
  dat <- simData(N, m, random_treatment)
  estimators <- value_estimators(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
                                 event_misspecification,
                                 sampling_misspecification,
                                 propensity_misspecification, 
                                 censor_misspecification,
                                 model_propensity,
                                 model_survival,
                                 model_calibration,
                                 crossfitV)
  # true value
  Vd_naive[i] <- Vd(estimators$d_naive)
  Vd_ipsw[i] <- Vd(estimators$d_ipsw)
  Vd_cwipw[i] <- Vd(estimators$d_cwipw)
  Vd_cwor[i] <- Vd(estimators$d_cwor)
  Vd_ort[i] <- Vd(estimators$d_ort)
  Vd_acw[i] <- Vd(estimators$d_acw)
  
  # proportion of correct decision
  pcd_naive[i] <- pcd(estimators$d_naive)
  pcd_ipsw[i] <- pcd(estimators$d_ipsw)
  pcd_cwipw[i] <- pcd(estimators$d_cwipw)
  pcd_cwor[i] <- pcd(estimators$d_cwor)
  pcd_ort[i] <- pcd(estimators$d_ort)
  pcd_acw[i] <- pcd(estimators$d_acw)
  
  Vhat_naive[i] <- estimators$Vnaive
  Vhat_ipsw[i] <- estimators$Vipsw
  Vhat_cwipw[i] <- estimators$Vcwipw
  Vhat_cwor[i] <- estimators$Vcwor
  Vhat_ort[i] <- estimators$Vort
  Vhat_acw[i] <- estimators$Vacw
  
  # bootstrap variance estimation
  # parallel computation
  # n.workers <- detectCores()
  cl <- makePSOCKcluster(48)
  clusterExport(cl, c("event_misspecification", "sampling_misspecification",
                      "propensity_misspecification", "censor_misspecification",
                      "dat",
                      "model_propensity", "model_survival",
                      "model_calibration", "crossfitV",
                      "simData", "compute_ps", "survivalFit",
                      "compute_RMST", "lamFun", "compute_calibration_weights",
                      "value_estimators", "Vd", "pcd", "L",
                      "largeXtarget", "largeT1L", "largeT0L", "Aopt"))
  registerDoParallel(cl)
  
  bootV <- foreach(b = 1:B, .combine = 'cbind',
                   .packages = c("rgenoud", "survival", "nleqslv")) %dopar% {
                     n <- nrow(dat$Xsource)
                     m <- nrow(dat$Xtarget)
                     idx_b <- sample(n, replace = TRUE)
                     estimator_b <- value_estimators(dat$Xsource[idx_b, ], dat$A[idx_b], dat$U[idx_b], dat$status[idx_b], 
                                                     dat$Xtarget[sample(m, replace = TRUE), ],
                                                     event_misspecification,
                                                     sampling_misspecification,
                                                     propensity_misspecification, 
                                                     censor_misspecification,
                                                     model_propensity,
                                                     model_survival,
                                                     model_calibration,
                                                     crossfitV)
                     c(estimator_b$Vnaive, estimator_b$Vipsw, estimator_b$Vcwipw,
                       estimator_b$Vcwor, estimator_b$Vort, estimator_b$Vacw)
                   }
  stopCluster(cl)
  # end
  sd.naive <- sd(bootV[1, ])
  sd.ipsw <- sd(bootV[2, ])
  sd.cwipw <- sd(bootV[3, ])
  sd.cwor <- sd(bootV[4, ])
  sd.ort <- sd(bootV[5, ])
  sd.acw <- sd(bootV[6, ])
  
  # results summary
  res.naive[i, ] <- c(Vhat_naive[i] - Vmax, sd.naive,
                      as.numeric(Vmax >= (Vhat_naive[i] - qnorm(0.975) * sd.naive)
                                 & Vmax <= (Vhat_naive[i] + qnorm(0.975) * sd.naive)),
                      as.numeric(Vmax >= quantile(bootV[1, ], 0.025)
                                 & Vmax <= quantile(bootV[1, ], 0.975)))
  res.ipsw[i, ] <- c(Vhat_ipsw[i] - Vmax, sd.ipsw,
                     as.numeric(Vmax >= (Vhat_ipsw[i] - qnorm(0.975) * sd.ipsw)
                                & Vmax <= (Vhat_ipsw[i] + qnorm(0.975) * sd.ipsw)),
                     as.numeric(Vmax >= quantile(bootV[2, ], 0.025)
                                & Vmax <= quantile(bootV[2, ], 0.975)))
  res.cw.ipw[i, ] <- c(Vhat_cwipw[i] - Vmax, sd.cwipw,
                       as.numeric(Vmax >= (Vhat_cwipw[i] - qnorm(0.975) * sd.cwipw)
                                  & Vmax <= (Vhat_cwipw[i] + qnorm(0.975) * sd.cwipw)),
                       as.numeric(Vmax >= quantile(bootV[3, ], 0.025)
                                  & Vmax <= quantile(bootV[3, ], 0.975)))
  res.cw.or[i, ] <- c(Vhat_cwor[i] - Vmax, sd.cwor,
                      as.numeric(Vmax >= (Vhat_cwor[i] - qnorm(0.975) * sd.cwor)
                                 & Vmax <= (Vhat_cwor[i] + qnorm(0.975) * sd.cwor)),
                      as.numeric(Vmax >= quantile(bootV[4, ], 0.025)
                                 & Vmax <= quantile(bootV[4, ], 0.975)))
  res.ort[i, ] <- c(Vhat_ort[i] - Vmax, sd.ort,
                    as.numeric(Vmax >= (Vhat_ort[i] - qnorm(0.975) * sd.ort)
                               & Vmax <= (Vhat_ort[i] + qnorm(0.975) * sd.ort)),
                    as.numeric(Vmax >= quantile(bootV[5, ], 0.025)
                               & Vmax <= quantile(bootV[5, ], 0.975)))
  res.acw[i, ] <- c(Vhat_acw[i] - Vmax, sd.acw,
                    as.numeric(Vmax >= (Vhat_acw[i] - qnorm(0.975) * sd.acw)
                               & Vmax <= (Vhat_acw[i] + qnorm(0.975) * sd.acw)),
                    as.numeric(Vmax >= quantile(bootV[6, ], 0.025)
                               & Vmax <= quantile(bootV[6, ], 0.975)))
  cat(i, as.character(Sys.time()), "\n")
}

result3 <- list(Vd = list(Vd_naive = Vd_naive,
                          Vd_ipsw = Vd_ipsw,
                          Vd_cwipw = Vd_cwipw,
                          Vd_cwor = Vd_cwor,
                          Vd_ort = Vd_ort,
                          Vd_acw = Vd_acw),
                pcd = list(pcd_naive = pcd_naive,
                           pcd_ipsw = pcd_ipsw,
                           pcd_cwipw = pcd_cwipw,
                           pcd_cwor = pcd_cwor,
                           pcd_ort = pcd_ort,
                           pcd_acw = pcd_acw),
                Vhat = list(Vhat_naive = Vhat_naive,
                            Vhat_ipsw = Vhat_ipsw,
                            Vhat_cwipw = Vhat_cwipw,
                            Vhat_cwor = Vhat_cwor,
                            Vhat_ort = Vhat_ort,
                            Vhat_acw = Vhat_acw),
                res = list(res.naive = res.naive,
                           res.ipsw = res.ipsw,
                           res.cw.ipw = res.cw.ipw,
                           res.cw.or = res.cw.or,
                           res.ort = res.ort,
                           res.acw = res.acw))

save(result3, file = "survITR3.RData")


#########
## CASE 4

Vd_naive <- Vd_ipsw <- Vd_cwipw <- Vd_cwor <- Vd_ort <- Vd_acw <- rep(0, n.MonteCarlo)
pcd_naive <- pcd_ipsw <- pcd_cwipw <- pcd_cwor <- pcd_ort <- pcd_acw <- rep(0, n.MonteCarlo)
Vhat_naive <- Vhat_ipsw <- Vhat_cwipw <- Vhat_cwor <- Vhat_ort <- Vhat_acw <- rep(0, n.MonteCarlo)

res.naive <- matrix(0, n.MonteCarlo, 4)
colnames(res.naive) <- c("bias", "se", "coverageW", "coverageQ")
res.ipsw <- res.cw.ipw <- res.cw.or <- res.ort <- res.acw <- res.naive

event_misspecification <- TRUE
sampling_misspecification <- TRUE
propensity_misspecification <- TRUE
censor_misspecification <- TRUE

set.seed(6)
for (i in 1:n.MonteCarlo) {
  dat <- simData(N, m, random_treatment)
  estimators <- value_estimators(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
                                 event_misspecification,
                                 sampling_misspecification,
                                 propensity_misspecification, 
                                 censor_misspecification,
                                 model_propensity,
                                 model_survival,
                                 model_calibration,
                                 crossfitV)
  # true value
  Vd_naive[i] <- Vd(estimators$d_naive)
  Vd_ipsw[i] <- Vd(estimators$d_ipsw)
  Vd_cwipw[i] <- Vd(estimators$d_cwipw)
  Vd_cwor[i] <- Vd(estimators$d_cwor)
  Vd_ort[i] <- Vd(estimators$d_ort)
  Vd_acw[i] <- Vd(estimators$d_acw)
  
  # proportion of correct decision
  pcd_naive[i] <- pcd(estimators$d_naive)
  pcd_ipsw[i] <- pcd(estimators$d_ipsw)
  pcd_cwipw[i] <- pcd(estimators$d_cwipw)
  pcd_cwor[i] <- pcd(estimators$d_cwor)
  pcd_ort[i] <- pcd(estimators$d_ort)
  pcd_acw[i] <- pcd(estimators$d_acw)
  
  Vhat_naive[i] <- estimators$Vnaive
  Vhat_ipsw[i] <- estimators$Vipsw
  Vhat_cwipw[i] <- estimators$Vcwipw
  Vhat_cwor[i] <- estimators$Vcwor
  Vhat_ort[i] <- estimators$Vort
  Vhat_acw[i] <- estimators$Vacw
  
  # bootstrap variance estimation
  # parallel computation
  # n.workers <- detectCores()
  cl <- makePSOCKcluster(48)
  clusterExport(cl, c("event_misspecification", "sampling_misspecification",
                      "propensity_misspecification", "censor_misspecification",
                      "dat",
                      "model_propensity", "model_survival",
                      "model_calibration", "crossfitV",
                      "simData", "compute_ps", "survivalFit",
                      "compute_RMST", "lamFun", "compute_calibration_weights",
                      "value_estimators", "Vd", "pcd", "L",
                      "largeXtarget", "largeT1L", "largeT0L", "Aopt"))
  registerDoParallel(cl)
  
  bootV <- foreach(b = 1:B, .combine = 'cbind',
                   .packages = c("rgenoud", "survival", "nleqslv")) %dopar% {
                     n <- nrow(dat$Xsource)
                     m <- nrow(dat$Xtarget)
                     idx_b <- sample(n, replace = TRUE)
                     estimator_b <- value_estimators(dat$Xsource[idx_b, ], dat$A[idx_b], dat$U[idx_b], dat$status[idx_b], 
                                                     dat$Xtarget[sample(m, replace = TRUE), ],
                                                     event_misspecification,
                                                     sampling_misspecification,
                                                     propensity_misspecification, 
                                                     censor_misspecification,
                                                     model_propensity,
                                                     model_survival,
                                                     model_calibration,
                                                     crossfitV)
                     c(estimator_b$Vnaive, estimator_b$Vipsw, estimator_b$Vcwipw,
                       estimator_b$Vcwor, estimator_b$Vort, estimator_b$Vacw)
                   }
  stopCluster(cl)
  # end
  sd.naive <- sd(bootV[1, ])
  sd.ipsw <- sd(bootV[2, ])
  sd.cwipw <- sd(bootV[3, ])
  sd.cwor <- sd(bootV[4, ])
  sd.ort <- sd(bootV[5, ])
  sd.acw <- sd(bootV[6, ])
  
  # results summary
  res.naive[i, ] <- c(Vhat_naive[i] - Vmax, sd.naive,
                      as.numeric(Vmax >= (Vhat_naive[i] - qnorm(0.975) * sd.naive)
                                 & Vmax <= (Vhat_naive[i] + qnorm(0.975) * sd.naive)),
                      as.numeric(Vmax >= quantile(bootV[1, ], 0.025)
                                 & Vmax <= quantile(bootV[1, ], 0.975)))
  res.ipsw[i, ] <- c(Vhat_ipsw[i] - Vmax, sd.ipsw,
                     as.numeric(Vmax >= (Vhat_ipsw[i] - qnorm(0.975) * sd.ipsw)
                                & Vmax <= (Vhat_ipsw[i] + qnorm(0.975) * sd.ipsw)),
                     as.numeric(Vmax >= quantile(bootV[2, ], 0.025)
                                & Vmax <= quantile(bootV[2, ], 0.975)))
  res.cw.ipw[i, ] <- c(Vhat_cwipw[i] - Vmax, sd.cwipw,
                       as.numeric(Vmax >= (Vhat_cwipw[i] - qnorm(0.975) * sd.cwipw)
                                  & Vmax <= (Vhat_cwipw[i] + qnorm(0.975) * sd.cwipw)),
                       as.numeric(Vmax >= quantile(bootV[3, ], 0.025)
                                  & Vmax <= quantile(bootV[3, ], 0.975)))
  res.cw.or[i, ] <- c(Vhat_cwor[i] - Vmax, sd.cwor,
                      as.numeric(Vmax >= (Vhat_cwor[i] - qnorm(0.975) * sd.cwor)
                                 & Vmax <= (Vhat_cwor[i] + qnorm(0.975) * sd.cwor)),
                      as.numeric(Vmax >= quantile(bootV[4, ], 0.025)
                                 & Vmax <= quantile(bootV[4, ], 0.975)))
  res.ort[i, ] <- c(Vhat_ort[i] - Vmax, sd.ort,
                    as.numeric(Vmax >= (Vhat_ort[i] - qnorm(0.975) * sd.ort)
                               & Vmax <= (Vhat_ort[i] + qnorm(0.975) * sd.ort)),
                    as.numeric(Vmax >= quantile(bootV[5, ], 0.025)
                               & Vmax <= quantile(bootV[5, ], 0.975)))
  res.acw[i, ] <- c(Vhat_acw[i] - Vmax, sd.acw,
                    as.numeric(Vmax >= (Vhat_acw[i] - qnorm(0.975) * sd.acw)
                               & Vmax <= (Vhat_acw[i] + qnorm(0.975) * sd.acw)),
                    as.numeric(Vmax >= quantile(bootV[6, ], 0.025)
                               & Vmax <= quantile(bootV[6, ], 0.975)))
  cat(i, as.character(Sys.time()), "\n")
}

result4 <- list(Vd = list(Vd_naive = Vd_naive,
                          Vd_ipsw = Vd_ipsw,
                          Vd_cwipw = Vd_cwipw,
                          Vd_cwor = Vd_cwor,
                          Vd_ort = Vd_ort,
                          Vd_acw = Vd_acw),
                pcd = list(pcd_naive = pcd_naive,
                           pcd_ipsw = pcd_ipsw,
                           pcd_cwipw = pcd_cwipw,
                           pcd_cwor = pcd_cwor,
                           pcd_ort = pcd_ort,
                           pcd_acw = pcd_acw),
                Vhat = list(Vhat_naive = Vhat_naive,
                            Vhat_ipsw = Vhat_ipsw,
                            Vhat_cwipw = Vhat_cwipw,
                            Vhat_cwor = Vhat_cwor,
                            Vhat_ort = Vhat_ort,
                            Vhat_acw = Vhat_acw),
                res = list(res.naive = res.naive,
                           res.ipsw = res.ipsw,
                           res.cw.ipw = res.cw.ipw,
                           res.cw.or = res.cw.or,
                           res.ort = res.ort,
                           res.acw = res.acw))

save(result4, file = "survITR4.RData")


############
## ML models

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
  lambdaC <- 0.2
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


#########################
## performance comparison

# simulation setup
N <- 6e5
m <- 24000
random_treatment <- FALSE

L <- 4

event_misspecification <- FALSE
sampling_misspecification <- FALSE
propensity_misspecification <- FALSE
censor_misspecification <- FALSE

model_propensity <- "RF"
model_survival <- "RF"
model_calibration <- "d2"
crossfitV <- 2

n.MonteCarlo <- 200

Vd_naive <- Vd_ipsw <- Vd_cwipw <- Vd_cwor <- Vd_ort <- Vd_acw <- rep(0, n.MonteCarlo)
pcd_naive <- pcd_ipsw <- pcd_cwipw <- pcd_cwor <- pcd_ort <- pcd_acw <- rep(0, n.MonteCarlo)
Vhat_naive <- Vhat_ipsw <- Vhat_cwipw <- Vhat_cwor <- Vhat_ort <- Vhat_acw <- rep(0, n.MonteCarlo)

set.seed(67)
for (i in 1:n.MonteCarlo) {
  dat <- simData(N, m, random_treatment)
  estimators <- value_estimators(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
                                 event_misspecification,
                                 sampling_misspecification,
                                 propensity_misspecification, 
                                 censor_misspecification,
                                 model_propensity,
                                 model_survival,
                                 model_calibration,
                                 crossfitV)
  # true value
  Vd_naive[i] <- Vd(estimators$d_naive)
  Vd_ipsw[i] <- Vd(estimators$d_ipsw)
  Vd_cwipw[i] <- Vd(estimators$d_cwipw)
  Vd_cwor[i] <- Vd(estimators$d_cwor)
  Vd_ort[i] <- Vd(estimators$d_ort)
  Vd_acw[i] <- Vd(estimators$d_acw)
  
  # proportion of correct decision
  pcd_naive[i] <- pcd(estimators$d_naive)
  pcd_ipsw[i] <- pcd(estimators$d_ipsw)
  pcd_cwipw[i] <- pcd(estimators$d_cwipw)
  pcd_cwor[i] <- pcd(estimators$d_cwor)
  pcd_ort[i] <- pcd(estimators$d_ort)
  pcd_acw[i] <- pcd(estimators$d_acw)
  
  Vhat_naive[i] <- estimators$Vnaive
  Vhat_ipsw[i] <- estimators$Vipsw
  Vhat_cwipw[i] <- estimators$Vcwipw
  Vhat_cwor[i] <- estimators$Vcwor
  Vhat_ort[i] <- estimators$Vort
  Vhat_acw[i] <- estimators$Vacw
  
  cat(i, as.character(Sys.time()), "\n")
}

save(Vd_naive, Vd_ipsw, Vd_cwipw, Vd_cwor, Vd_ort, Vd_acw,
     pcd_naive, pcd_ipsw, pcd_cwipw, pcd_cwor, pcd_ort, pcd_acw,
     Vhat_naive, Vhat_ipsw, Vhat_cwipw, Vhat_cwor, Vhat_ort, Vhat_acw,
     Vmax,
     file = "survITR_NP1.RData")


##############
## sample size

# proposed estimator
ACW <- function(Xsource, A, U, status, Xtarget,
                event_misspecification = FALSE,
                sampling_misspecification = FALSE,
                propensity_misspecification = FALSE, 
                censor_misspecification = FALSE,
                model_propensity = "RF",
                model_survival = "RF",
                model_calibration = "d2",
                crossfitV = 2) {
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
  
  # ACW
  RMSTort0 <- compute_RMST(sF$S0_new, U, L)
  RMSTort1 <- compute_RMST(sF$S1_new, U, L)
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
  
  dA1 <- as.numeric(cbind(1, Xtarget) %*% d_acw > 0)
  dA2 <- as.numeric(cbind(1, Xsource) %*% d_acw > 0)
  SEacw <- sd(dA1 * RMSTort1 + (1 - dA1) * RMSTort0) / sqrt(m) + sd(as.numeric(A == dA2) * q * RMSTaug) * sqrt(n)
  CPacw <- as.numeric(Vmax > (Vacw - qnorm(0.975) * SEacw)
                      & Vmax < (Vacw + qnorm(0.975) * SEacw))
  return(list(Vacw = Vacw,
              CPacw = CPacw,
              SEacw = SEacw))
}

n.MonteCarlo <- 200
Vacw1 <- Vacw2 <- Vacw3 <- Vacw4 <- Vacw5 <- Vacw6 <- rep(NA, n.MonteCarlo)
CPacw1 <- CPacw2 <- CPacw3 <- CPacw4 <- CPacw5 <- CPacw6 <- rep(NA, n.MonteCarlo)
SEacw1 <- SEacw2 <- SEacw3 <- SEacw4 <- SEacw5 <- SEacw6 <- rep(NA, n.MonteCarlo)

set.seed(66)
for (i in 1:n.MonteCarlo) {
  dat <- simData(5e4, 2000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw1[i] <- acw$Vacw
  CPacw1[i] <- acw$CPacw
  SEacw1[i] <- acw$SEacw
  
  dat <- simData(1e5, 4000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw2[i] <- acw$Vacw
  CPacw2[i] <- acw$CPacw
  SEacw2[i] <- acw$SEacw
  
  dat <- simData(2e5, 8000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw3[i] <- acw$Vacw
  CPacw3[i] <- acw$CPacw
  SEacw3[i] <- acw$SEacw
  
  dat <- simData(4e5, 16000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw4[i] <- acw$Vacw
  CPacw4[i] <- acw$CPacw
  SEacw4[i] <- acw$SEacw
  
  dat <- simData(6e5, 24000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw5[i] <- acw$Vacw
  CPacw5[i] <- acw$CPacw
  SEacw5[i] <- acw$SEacw
  
  dat <- simData(8e5, 32000, FALSE)
  acw <- ACW(dat$Xsource, dat$A, dat$U, dat$status, dat$Xtarget,
             FALSE, FALSE, FALSE, FALSE,
             "RF", "RF", "d2", 2)
  Vacw6[i] <- acw$Vacw
  CPacw6[i] <- acw$CPacw
  SEacw6[i] <- acw$SEacw
  
  cat(i, as.character(Sys.time()), "\n")
}

save(Vacw1, Vacw2, Vacw3, Vacw4, Vacw5, Vacw6,
     CPacw1, CPacw2, CPacw3, CPacw4, CPacw5, CPacw6,
     SEacw1, SEacw2, SEacw3, SEacw4, SEacw5, SEacw6,
     Vmax,
     file = "survITR_NP2.RData")

