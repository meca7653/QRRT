#######################################
#' @name QRRT
#' @aliases QRRT
#' @title QRRT
#' @description Analysis QRRT data
#'
#' @param Formula same model specification as linear model (lm)
#' @param Data a data frame containing the variables in the model
#' @param Disperse standard deviation of the random starting values for beta
#' @param beta specified starting value for regression coefficients
#' @param n_times number of random starting values for beta
#' @param b_distribution prob distribution on 0, 1, ..., m, m + 1
#' @import stats
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
QRRT <- function(Formula, Data, Disperse = 1, beta = NULL, n_times = 1,
                offset = rep(0, nobs),
                b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50){
  m <- length(b_distribution) - 2 # b distribution has mass on 0, 1, ..., m and m+1
  #--------------------------------------------------------------------------------
  # Conditional expectation for the E-step.
  expectation <- function(x, Ri, beta, m, b_distribution, offset){
    lambda <- exp(x %*% beta + offset)
    expectation <- rep(NA, length(Ri))
    for (i in c(1:length(Ri))){
      if (Ri[i] < (m + 1)){
        expectation[i] <- dpois(Ri[i], lambda[i]) * b_distribution[m + 2] /
          (b_distribution[Ri[i] + 1] + dpois(Ri[i], lambda[i]) * b_distribution[m + 2])
      }else{
        expectation[i] <- 1
      }
    }
    return(expectation)
  }
  #--------------------------------------------------------------------------------
  ##  Likelihood function
  Likelihood <- function(x, Ri, beta, m, b_distribution, offset){
    ex <- expectation(x, Ri, beta, m,b_distribution, offset)
    lambda <- exp(x %*% beta + offset)
    likelihood <- -sum(ex * lambda) + sum(x %*% beta * Ri * ex)
    return(likelihood)
  }
  #--------------------------------------------------------------------------------
  N <- dim(Data)[1]
  a <- model.frame(Formula, Data)
  Y <- model.response(a)
  x.matrix <- model.matrix(Formula, data = a)
  Data <- a
  x_max <- apply(abs(x.matrix), 2, max)
  if(sum(x_max == 0) >0){
    print(paste0("Error: singularity in design matrix: column ", which(x_max == 0),
                 " identically zero"))
    return()
  }
  x <- x.matrix %*% (diag(1 / x_max))

  n <- dim(x)[1]

  #######################
  Ri <- model.response(a)
  l1_result <- c()

  beta_initial <- matrix(NA, nrow = dim(x)[2], ncol = n_times)

  beta_1_result <- matrix(NA, nrow = dim(x)[2], ncol = n_times)
  if(is.null(beta)){
    beta <- rnorm(n = dim(x)[2], sd = Disperse)# random starting values
  }else{
    beta <- beta
  }

  for (iii in c(1:n_times)){
    beta_initial[, iii] <- beta
    #-----------------------------------------------------------------------------
    # Initialize the EM algorithm, starting from random beta.
    beta <- rnorm(n = length(beta), beta, sd = Disperse)
    l1 <- Likelihood(x, Ri, beta, m, b_distribution, offset)
    weight <- expectation(x, Ri, beta, m, b_distribution, offset)
    data1 <- data.frame(Data, weight)
    lm1 <- glm.fit(x, Y, family = poisson(), weights = weight, offset = offset)
    beta_1 <- lm1$coefficients
    if (sum(is.na(beta_1))>0){
      beta_1[which(is.na(beta_1))] = 0
    }
    l2 <- Likelihood(x, Ri, beta_1, m, b_distribution, offset)
    #-----------------------------------------------------------------------------
    # Iterate EM algorithm to convergence.
    j <- 1
    while(abs(l1 - l2) > 1e-8){
      l1 <- l2
      # print(l1)
      beta <- beta_1
      #-----------------------------------------------------------------------------
      # E-step of EM algorithm: compute weights under current maximized model.
      weight <- expectation(x, Ri, beta, m, b_distribution, offset)
      # M-step of EM algorithm:  maximize the weighted log-likelihood
      # for Poisson regression, using case weights obtained from the E-step.
      lm1 <- glm.fit(x, Y, family = poisson(), weights = weight, offset = offset)
      beta_1 <- lm1$coefficients
      if (sum(is.na(beta_1)) > 0){
        beta_1[which(is.na(beta_1))] <- 0
      }
      j <- j + 1
      # Evaluate likelihood of current estimates, to assess convergence.
      l2 <- Likelihood(x, Ri, beta_1, m, b_distribution, offset)
    }
    # End of EM algorithm for this random start of EM algorithm.
    #-----------------------------------------------------------------------------
    # Store log-likelihoods for each raandom start, for later comparison.
    l1_result <- c(l1_result, l1)
    beta_1_result[,iii] <- beta_1
    print(iii)
    beta <- beta_1
  }
  # End of loop across n_times random starts of EM algorithm.
  #-----------------------------------------------------------------------------
  beta_name <- names(beta_1)
  beta_1 <- beta_1_result[, which.max(l1_result)]
  names(beta_1) <- beta_name
  lambda <- exp(x %*% beta_1 + offset)
  fish <- rep(NA, length(Ri))
  # fr <- rep(NA, length(Ri))
  # dfr <- matrix(NA, nrow = length(Ri), ncol = length(beta_1))
  R <- (0:m)
  # ###############
  # for (i in (1:length(Ri))){
  #   if (sum(R == Ri[i]) > 0){
  #     fr[i] <- b_distribution[m + 2] * dpois(Ri[i], lambda[i]) +
  #       b_distribution[which(R == Ri[i])]
  #   }else{
  #     fr[i] = b_distribution[m + 2] * dpois(Ri[i], lambda[i])
  #   }
  #
  #   for (j in c(1:length(beta_1))){
  #     dfr[i, j] <- x[i, j] * 1 / fr[i] * b_distribution[m + 2] *
  #       (Ri[i] - lambda[i]) * dpois(Ri[i], lambda[i])
  #   }
  #
  # }
  # dfr_final <- data.frame(colSums(dfr))
  for (i in (1:length(Ri))){
    fr1 <- b_distribution[m + 2] * dpois(R, lambda[i]) + b_distribution[1:(m + 1)]
    fish[i] <- b_distribution[m + 2] *
      (-lambda[i] + sum(dpois(R, lambda[i]) * b_distribution[1:(m + 1)] *
                          (R - lambda[i]) ^ 2 / fr1))
  }
  fish_info <- matrix(rep(NA, length(beta_1) ^ 2), ncol = length(beta_1))
  for (i in (1:length(beta_1))){
    for (j in (1:length(beta_1))){
      fish_info[i, j] <- -sum(fish * x[, i] * x[, j])
    }
  }
  inverse_fisher_info <- solve(fish_info)
  f <- function(x, beta, Ri){
    #function for likelihood
    lambda <- exp(x%*%beta)
    f_ri <- rep(NA, length(Ri))
    for (i in c(1:length(Ri))){
      if (Ri[i] <= m){
        f_ri[i] <- b_distribution[Ri[i]+1] + b_distribution[m+2]*dpois( x = Ri[i], lambda = lambda[i])
      }else{
        f_ri[i] <- b_distribution[m+2]*dpois(x = Ri[i], lambda = lambda[i])
      }


    }
    return(f_ri)

  }


  fr <- F_Ri <- f(x, beta, Ri)

  first_div <- function(x, beta, Ri, f_ri){
    #first devirative for log-likelihood
    lambda <- exp(x%*%beta)
    first_div <- matrix(rep(NA, length(Ri)*length(beta)), ncol = length(beta))
    for (i in c(1:length(Ri))){
      for (j in c(1:length(beta))){
        first_div[i, j] <- 1/f_ri[i] *x[i,j]*b_distribution[m+2]*(-dpois(x = Ri[i], lambda = lambda[i])*lambda[i] + Ri[i]*dpois(x = Ri[i], lambda = lambda[i]))

      }

    }
    return(first_div)
  }

  fir_div <- (first_div(x = x, beta = beta, Ri = Ri, f_ri = F_Ri))
  dfr_final <- data.frame(colSums(fir_div))

  second_div <- function(x, beta, Ri, f_ri, fir_div){
    #second derivative for log-likelihood
    second <- matrix(rep(NA, length(beta)^2), ncol = length(beta))
    lambda <- exp(x%*%beta)
    for (j in c(1:length(beta))){
      for (k in c(1:length(beta))){
        second[j,k] <- sum(-fir_div[,j]*fir_div[,k] + 1/f_ri*x[,j]*x[,k]*b_distribution[m+2]*(dpois(x = Ri, lambda))*(lambda^2 - (Ri+1)*lambda + Ri^2 - Ri*lambda))
      }
    }
    return(second)
  }
  sec_div <- second_div(x, beta = beta, Ri = Ri, f_ri = F_Ri,  fir_div = fir_div)

  fisher_est <- diag(1 / x_max) %*% solve(sec_div) %*% (diag(1 / x_max))



  p_value <- 2 * pnorm(-abs(beta_1) / sqrt(diag(inverse_fisher_info)))
  out_pre <- cbind(beta_1, sqrt(diag(inverse_fisher_info)), beta_1 / sqrt(diag(inverse_fisher_info)), p_value)
  out <- cbind(beta_1 / x_max, sqrt(diag(inverse_fisher_info)) / x_max,
               beta_1 / sqrt(diag(inverse_fisher_info)), p_value)
  out <- data.frame(out)
  names(out) <- c("Estimate", "Std.Error", "t-statistic", "Pr(>|t|)")
  max_like <- sum(log(fr))
  AIC  <- -2 * max_like + 2 * length(beta_1)
  names(AIC) <- "AIC"
  missing <- N - n
  sample_size <- n
  #--------------------------------------------------------------------------------
  #######Wald test########
  beta_length <- length(beta_1)
  if (beta_length >1){
    W <- t(beta_1[2:beta_length]) %*% fish_info[2:beta_length, 2:beta_length] %*%
      beta_1[2:beta_length]
    p_value <- 1 - pchisq(W, df = (beta_length - 1))
  }else{
    W <- t(beta_1) %*% fish_info %*% beta_1
    p_value <- 1 - pchisq(W, df = (beta_length))
  } # Wald test does not seem to be returned.  Might delete.
  #--------------------------------------------------------------------------------
  ########output#####

  inverse_fisher_info_report <- diag(1 / x_max) %*% inverse_fisher_info %*% (diag(1 / x_max))
  result <- list(out, AIC,max_like, missing, n, l1_result, dfr_final,
                 inverse_fisher_info_report, out_pre, fisher_est)
  names(result) <- c("out", "AIC", "Maximized_Log_Likelihood",
                     "Number_Missing_Values",
                     "Sample_Size", "l1_result", "score", "covariance", "out_pre", "fisher_est")
  return(result)
}
#####End QRRT function.
#--------------------------------------------------------------------------------



