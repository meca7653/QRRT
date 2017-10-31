rm(list = ls())
library(foreach)
source("R/simulation_QRRT.R")
beta <- cbind(c(1.5, 1.0, -0.5, 0.4, 0.3, 0.2))
n <- 400 # sample size
set.seed(12345)
x1 <- sample(15:70, size = n, replace = TRUE) / 50
x2 <- rnorm(n, 0, 1)
x3 <- sample(c("A", "B", "C"), size = n, replace = TRUE)
dummyB <- 1 * (x3 == "B")
dummyC <- 1 * (x3 == "C")
X <- cbind(rep(1, n), x1, x2, dummyB, dummyC, x1 * x2)
# Using same known distribution as Conteh et al. (2015).
# This is a distribution on 0, 1, ..., m, m + 1.
b_distribution <- c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50
sum(b_distribution)
m <- length(b_distribution) - 2

#-----------------simulation function for 1000 times (no need to run)--------------------
# res_full = foreach(ii = c(1:1000), .combine = "rbind")%dopar%{
#
#   lambda <- exp( X %*% beta)
#   Ti <- rpois(n, lambda) # true response T_i
#   range(Ti)
#
#   # Draw the random integers B_i ~ b(r) for i = 1, 2, ..., n.
#   Bi <- sample(0:(m + 1), size = n, replace = TRUE, prob = b_distribution)
#   # If B_i < m + 1, observe B_i; otherwise, observe truth.
#   Ri <- Bi * (Bi < (m + 1)) + Ti * (Bi == (m + 1))
#   Sim_Data <- data.frame(x1, x2, x3, Ri)
#
#   fit_2way <- QRRT(Formula = Ri ~ (x1 + x2 + as.factor(x3)) ^ 2,
#                    Data = Sim_Data, Disperse = 1, beta = NULL, n_times = 1,
#                    offset = rep(0, n),
#                    b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50)
#   fit_2way$out
#   #------------------------------------------------------------------------------------
#   # Fit model correctly-specified as the actual simulation model.
#   fit_truemodel <- QRRT(Formula = Ri ~ (x1 + x2) ^ 2 + as.factor(x3),
#                         Data = Sim_Data, Disperse = 1, beta = NULL, n_times = 1,
#                         offset = rep(0, n),
#                         b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50)
#   fit_truemodel$out
#
#
#   list(Sim_Data, fit_2way , fit_truemodel)
#
#
# }
# save(res_full, file = "res_full_10_31.rda")


#-------------------load saved simulation result------------
load("res_full_10_31.rda")
#-------------fisher calculation-----------------------
lambda <- exp(X%*%beta)
R <- (0:m)
fish <- rep(NA, n)
for (i in (1:n)){
  fr1 <- b_distribution[m+2]*dpois(R,lambda[i]) + b_distribution[1:(m+1)]
  fish[i] <- b_distribution[m+2]*(-lambda[i] + sum(dpois(R, lambda[i])*b_distribution[1:(m+1)]*(R - lambda[i])^2/fr1))
}
fish_info <- matrix(rep(NA, length(beta)^2), ncol = length(beta))
for (i in (1:length(beta))){
  for (j in (1:length(beta))){
    fish_info[i,j] <- -sum(fish*X[,i]*X[,j])
  }
}
fisher <- solve(fish_info)
diag(fisher)

#--------------analysis simulation result------------------
head(res_full)

Sim_Data <- res_full[,1][[1]]
a <- model.frame(Ri ~ (x1 + x2) ^ 2 + as.factor(x3), Sim_Data)
Y <- model.response(a)
x.matrix <- model.matrix(Ri ~ (x1 + x2) ^ 2 + as.factor(x3), data = a)

a_new <- model.frame(Ri ~ (x1 + x2 + as.factor(x3)) ^ 2, Sim_Data)
x.matrix <- model.matrix(Ri ~ (x1 + x2 + as.factor(x3)) ^ 2, data = a_new)

res_over_m = do.call(rbind, res_full[,2])
res_true_m = do.call(rbind, res_full[,3])

fisher_res_over_m = do.call(rbind, lapply(c(1:1000),
                                   FUN = function(ii) {diag((res_over_m[,8])[[ii]])}))
fisher_res_true_m = do.call(rbind, lapply(c(1:1000),
                                   FUN = function(ii) {diag((res_true_m[,8])[[ii]])}))

fisher_est_over_m = do.call(rbind, lapply(c(1:1000),
                                          FUN = function(ii) {-diag((res_over_m[,10])[[ii]])}))
fisher_est_true_m = do.call(rbind, lapply(c(1:1000),
                                          FUN = function(ii) {-diag((res_true_m[,10])[[ii]])}))

est_res_over_m = do.call(rbind, lapply(c(1:1000),
                                       FUN = function(ii) {((res_over_m[,1])[[ii]])[,1]}))

est_res_true_m = do.call(rbind, lapply(c(1:1000),
                                       FUN = function(ii) {((res_true_m[,1])[[ii]])[,1]}))

#-------------True model result---------------------
apply(est_res_true_m, 2, mean) # coef
apply(est_res_true_m, 2, var) # empirical variance
apply(fisher_res_true_m, 2, mean) # estimated variance using fisher
apply(fisher_est_true_m, 2, mean) # estimated variance using negative inverse second derivative
diag(fisher) # true fisher


#------------------ oversized model results ---------------------
apply(est_res_over_m, 2, mean) # coef
apply(est_res_over_m, 2, var) # empirical variance
apply(fisher_res_over_m, 2, mean) # estimated fisher
apply(fisher_est_over_m, 2, mean) # estimated variance using negative inverse second derivative


#----------------likelihood ratio test --------------
fisher_like_over_m = do.call(rbind, lapply(c(1:1000),
                                          FUN = function(ii) {((res_over_m[,3])[[ii]])}))
fisher_like_true_m = do.call(rbind, lapply(c(1:1000),
                                          FUN = function(ii) {((res_true_m[,3])[[ii]])}))
p_value = 1 - pchisq(- 2 * fisher_like_true_m + 2 * fisher_like_over_m, df = 4)
par(mfrow = c(1,2))
hist(p_value)
qqplot(qchisq(ppoints(100), df = 4), - 2 * fisher_like_true_m + 2 * fisher_like_over_m)
abline(0,1, col = "red")

sum(p_value < 0.05)
