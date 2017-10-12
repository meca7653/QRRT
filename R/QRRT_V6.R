#######################################
#' @name QRRT
#' @aliases QRRT
#' @title QRRT
#' @description Analysis QRRT data
#'
#' @param Formula same format as linear model (lm)
#' @param Data need a data frame
#' @param Disperse standard deviation of the random starting value
#' @param beta starting value
#' @param n_times number of random starting values
#' @param m max number on green balls
#' @param n_orange number of orange balls
#' @param n_green number of green balls
#' @import MASS
#' @import stats
#' @import mvtnorm
#'
#'
#' @return Data frame contains estimated beta, hypothesis result
#'
#' @examples
#' true_beta = c(log(5), 0, log(4), log(2), log(1))
#' n = 250 #sample size
#' age = 15:70
#' Age = sample(age, size = n, replace = TRUE)/100
#' continuous_x = rnorm(n, 0, 1)
#' discrete_x = sample(c("A", "B", "C"), size = n, replace = TRUE)
#' x = matrix(c(rep(1,n), Age, continuous_x, Age*continuous_x, as.factor(discrete_x)), ncol = 5)
#' Data = simulation_fun(n = 250, m = 8, n_orange = 25,n_green = c(6,7,4,2,2,1,1,1,1),
#' x = x, true_beta = true_beta)
#' res = QRRT(Formula = Y~.-1, Data = Data, Disperse = 1, beta = NULL, n_times = 1)
#'
#' @author Meng Cao, F. Jay Breidt
#' @export QRRT
QRRT = function(Formula, Data, Disperse = 1, beta = NULL, n_times = 1,
                m = 8, n_orange = 25,
                n_green = c(6,7,4,2,2,1,1,1,1)){
  ############################################
  expectation <- function(x, Ri, beta, m, p_balls){
    lambda = exp(x%*%beta)
    expectation = rep(NA, length(Ri))
    for (i in c(1:length(Ri))){
      if (Ri[i] <(m+1)){
        expectation[i] = dpois(Ri[i], lambda[i])*p_balls[m+2]/(p_balls[Ri[i]+1] + dpois(Ri[i], lambda[i])*p_balls[m+2])
      }else{
        expectation[i] = 1
      }
    }
    return(expectation)
  }
  ############Likelihood Function############
  Likelihood = function(x, Ri, beta, m, p_balls){
    ex = expectation(x, Ri, beta, m,p_balls)
    lambda = exp(x%*%beta)
    likelihood = -sum(ex*lambda) + sum(x%*%beta*Ri*ex)
    return(likelihood)
  }

  N<-dim(Data)[1]
  a <- model.frame(Formula,Data)
  Y<-model.response(a)
  x.matrix <- model.matrix(Formula, data=a)
  Data<-a
  x_max <- apply(abs(x.matrix), 2, max)
  if(sum(x_max == 0) >0){
    print(paste0("Error: singularity in design matrix: column ", which(x_max == 0),
                 " identically zero"))
    return()
  }
  x <- x.matrix %*% (diag(1/x_max))

  n <- dim(x)[1]

  ######################### balls ##########################
  balls = 0:(m+1)
  n_balls = sum(n_green) + n_orange
  p_balls = c(n_green, n_orange)/n_balls
  #######################
  Ri = model.response(a)
  l1_result = c()

  beta_initial <- matrix(NA, nrow = dim(x)[2], ncol = n_times)

  beta_1_result <- matrix(NA, nrow = dim(x)[2], ncol = n_times)
  if(is.null(beta)){
    beta<- rnorm(n = dim(x)[2], sd = Disperse)# random starting values
  }else{
    beta = beta
  }

  for (iii in c(1:n_times)){
    beta_initial[,iii] = beta
    #################normalization##############
    beta <- rnorm(n = length(beta),beta, sd = Disperse)
    l1 = Likelihood(x,Ri, beta, m, p_balls)
    weight = expectation(x, Ri, beta, m, p_balls)
    data1 = data.frame(Data, weight)
    lm1 = glm.fit(x,Y, family = poisson(), weights = weight)
    beta_1 = lm1$coefficients
    if (sum(is.na(beta_1))>0){
      beta_1[which(is.na(beta_1))] = 0
    }
    l2 = Likelihood(x,Ri, beta_1, m, p_balls)
    j = 1
    while(abs(l1 - l2) > 1e-8){
      l1 = l2
      # print(l1)
      beta = beta_1
      # beta_result = c(beta_result, beta[1])
      weight = expectation(x, Ri, beta, m, p_balls)
      lm1 = glm.fit(x,Y, family = poisson(), weights = weight)
      beta_1 = lm1$coefficients
      if (sum(is.na(beta_1))>0){
        beta_1[which(is.na(beta_1))] = 0
      }
      j = j+1
      l2 = Likelihood(x,Ri, beta_1, m, p_balls)

    }
    l1_result = c(l1_result, l1)
    beta_1_result[,iii] = beta_1
    print(iii)
    beta = beta_1
  }
  beta_name = names(beta_1)
  beta_1 = beta_1_result[,which.max(l1_result)]
  names(beta_1) = beta_name
  lambda = exp(x%*%beta_1)
  fish = rep(NA, length(Ri))
  fr = rep(NA, length(Ri))

  dfr = matrix(NA, nrow = length(Ri), ncol = length(beta_1))

  R = (0:m)
  ###############
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

  p_value<-2*pnorm(-abs(beta_1)/sqrt(diag(fisher)))
  out_pre <- cbind(beta_1 ,sqrt(diag(fisher)) , beta_1/sqrt(diag(fisher)), p_value)
  out = cbind(beta_1 /x_max, sqrt(diag(fisher)) / x_max, beta_1/sqrt(diag(fisher)), p_value)
  out = data.frame(out)
  names(out) = c("Estimate", "Std.Error", "t value" ,"Pr(>|t|)")
  max_like = sum(log(fr))


  AIC  = -2*max_like+2*length(beta_1)
  names(AIC) = "AIC"
  missing<-N-n
  sample_size<-n




  #######Wald test########

  beta_length = length(beta_1)
  if (beta_length >1){
    W = t(beta_1[2:beta_length])%*%fish_info[2:beta_length, 2:beta_length]%*%beta_1[2:beta_length]
    p_value = 1 - pchisq(W, df = (beta_length -1))
  }else{
    W = t(beta_1)%*%fish_info%*%beta_1
    p_value = 1 - pchisq(W, df = (beta_length))
  }


  ########output#####

  fisher_report <- diag(1/x_max) %*% fisher %*% (diag(1/x_max))
  result = list(out, AIC,max_like,missing,n, l1_result, dfr_final, fisher_report, out_pre)
  names(result) = c("out", "AIC","Maximized_Log_Likelihood","Number_Missing_Values","Sample_Size", "l1_result", "score","covariance", "out_pre")
  result
  return(result)
}

#####End QRRT function.



