library(foreach)
true_beta = c(log(5), log(3), log(4), log(2), log(5), 0.6)
n = 250 #sample size
continuous_x1 = sample(15:70, size = n, replace = TRUE)/50
continuous_x2 = rnorm(n, 0, 1)
discrete_x3 = sample(c("A", "B", "C"), size = n, replace = TRUE)
x = cbind(rep(1,n), continuous_x1, continuous_x2,
          continuous_x1*continuous_x2, discrete_x3 == "B", discrete_x3 == "C")

# colnames(x) = c("Intercept", "continuous_x1", "continuous_x2", "X1.X2", "discrete_x")
Data = simulation_fun(n = 250, m = 8, n_orange = 25,n_green = c(6,7,4,2,2,1,1,1,1),
                      x = x, true_beta = true_beta, offset = rep(3, n))

# names(Data) <- c("Y", paste0(rep("X", 5), c(1:5)))
# Data <- data.frame(Data$Y, continuous_x1, continuous_x2, discrete_x3)
names(Data)[1] = "R"
res = QRRT(Formula = R ~.-1
           , Data = Data, Disperse = 1, beta = NULL, n_times = 1, offset = rep(3, n))
res_full = foreach(ii = c(1:1000), .combine = "rbind")%dopar%{
  Data = simulation_fun(n = 250, m = 8, n_orange = 25,n_green = c(6,7,4,2,2,1,1,1,1),
                        x = x, true_beta = true_beta, offset = rep(3, n))
  # Data <- data.frame(Data$Y, continuous_x1, continuous_x2, discrete_x3)
  names(Data)[1] = "R"
  res = QRRT(Formula = R ~ . -1,
             Data = Data, Disperse = 1, beta = NULL, n_times = 1, offset = rep(3, n))

  Ri = Data$R
  offset = rep(3, n)
  lambda = exp(x%*%true_beta + offset)
  fish = rep(NA, length(Ri))
  fr = rep(NA, length(Ri))
  m = 8
  R = (0:m)
  balls = 0:(m+1)
  n_orange = 25
  n_green = c(6,7,4,2,2,1,1,1,1)
  n_balls = sum(n_green) + n_orange
  p_balls = c(n_green, n_orange)/n_balls
  ###############
  fish = rep(NA, length(Ri))
  fr = rep(NA, length(Ri))
  beta_1 = true_beta
  dfr = matrix(NA, nrow = length(Ri), ncol = length(beta_1))
  for (i in (1:length(Ri))){
    if (sum(R == Ri[i]) >0){
      fr[i] = p_balls[m+2]*dpois(Ri[i], lambda[i]) + p_balls[which(R == Ri[i])]
    }else{
      fr[i] = p_balls[m+2]*dpois(Ri[i], lambda[i])
    }

    for (j in c(1:length(beta_1))){
      dfr[i,j] = x[i,j]*1/fr[i]*p_balls[m+2]*(Ri[i] - lambda[i])*dpois(Ri[i], lambda[i])
    }

  }
  dfr_final = data.frame(colSums(dfr))
  for (i in (1:length(Ri))){
    fr1 = p_balls[m+2]*dpois(R,lambda[i]) + p_balls[1:(m+1)]
    fish[i] = p_balls[m+2]*(-lambda[i] + sum(dpois(R, lambda[i])*p_balls[1:(m+1)]*(R - lambda[i])^2/fr1))
  }
  fish_info = matrix(rep(NA, length(beta_1)^2), ncol = length(beta_1))
  for (i in (1:length(beta_1))){
    for (j in (1:length(beta_1))){
      fish_info[i,j] = -sum(fish*x[,i]*x[,j])
    }
  }
  fisher = solve(fish_info)
  list(res = res, fisher = fisher)


}

res_res <- do.call(rbind, res_full[, 1])
res_out <- do.call(rbind, res_res[,1])
save.image("res_full.rda")
apply(matrix(res_out[,1], nrow = 6), 1, mean)
apply(matrix(res_out[,1], nrow = 6), 1, var)
