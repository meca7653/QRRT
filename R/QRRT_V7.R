#######################################
#' @name QRRT
#' @aliases QRRT
#' @title Poisson regression for QRRT data
#' @description Poisson regression for quantitative randomized response technique (QRRT) data, based on maximum likelihood implemented via the EM algorithm
#'
#' @param Formula symbolic description of the model to be fitted using the same
#' specification as the lm function for linear models
#' @param Data data frame containing the variables in the model
#' @param beta NULL for random starting values or a specified starting value
#' for regression coefficients
#' @param Disperse standard deviation of the random starting values for beta
#' @param n_times number of random starting values for beta
#' @param offset NULL or a numeric vector of length equal to the number of cases
#' @param b_distribution probability distribution on 0, 1, ..., m, m + 1.
#' Default is from Conteh et al. (2015).
#' @param tolerance EM convergence criterion for change in log-likelihood
#' @import stats
#'
#'
#' @return List containing the following elements:
#' @param Results a matrix including maximum likelihood estimates of regression coefficients and  estimated asymptotic standard errors,
#' along with corresponding t-statistics and p-values
#'
#' @param AIC Akaike's An Information Criterion:
#' minus twice the maximized log-likelihood plus twice the number of parameters
#'
#' @param Maximized_Log_Likelihood log-likelihood at the estimated coefficient values
#' @param Likelihood_by_Starting_Value maximized likelihood at each of the n_times starting values
#' @param Number_Missing_Values number of missing value in Data
#' @param Score score vector evaluated at the maximum likelihood estimate of regression coefficient vector
#'
#' @param Inverse_Fisher_Information estimated variance-covariance matrix for maximum likelihood estimate of regression coefficient vector
#'
#'
#' @examples
#' #-------------------------------------------------------------------------------------
#' # Simulate a QRRT data set:
#' #-------------------------------------------------------------------------------------
#' # True count is Poisson with mean lambda =
#' # exp(beta_0+beta_1*x1+beta_2*x2+beta_3*dummyB+beta_4*dummyC+beta_5*x1*x2)
#' # where we have continuous predictors x1 and x2 and a categorical predictor
#' # x3 with levels A, B, C.  The categorical variable is encoded
#' # as a dummy variable for level B and a dummy variable for level C,
#' # with level A as baseline.
#' beta <- cbind(c(1.5, 1.0, -0.5, 0.4, 0.3, 0.2))
#' n <- 400 # sample size
#' set.seed(12345)
#' x1 <- sample(15:70, size = n, replace = TRUE) / 50
#' x2 <- rnorm(n, 0, 1)
#' x3 <- sample(c("A", "B", "C"), size = n, replace = TRUE)
#' dummyB <- 1 * (x3 == "B")
#' dummyC <- 1 * (x3 == "C")
#' X <- cbind(rep(1, n), x1, x2, dummyB, dummyC, x1 * x2)
#' lambda <- exp( X %*% beta)
#' Ti <- rpois(n, lambda) # true response T_i
#' range(Ti)
#' # Using same known distribution as Conteh et al. (2015).
#' # This is a distribution on 0, 1, ..., m, m + 1.
#' b_distribution <- c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50
#' sum(b_distribution)
#' m <- length(b_distribution) - 2
#' # Draw the random integers B_i ~ b(r) for i = 1, 2, ..., n.
#' Bi <- sample(0:(m + 1), size = n, replace = TRUE, prob = b_distribution)
#' # If B_i < m + 1, observe B_i; otherwise, observe truth.
#' Ri <- Bi * (Bi < (m + 1)) + Ti * (Bi == (m + 1))
#' Sim_Data <- data.frame(x1, x2, x3, Ri)
#' # Simulated data are now complete.
#'
#' #------------------------------------------------------------------------------------
#' # Now use the QRRT code to fit models to the simulated data.
#' #------------------------------------------------------------------------------------
#' # Fit model with all two-way interactions, using 10 random starting values.
#' fit_2way <- QRRT(Formula = Ri ~ (x1 + x2 + as.factor(x3)) ^ 2,
#'                  Data = Sim_Data, Disperse = 1, beta = NULL, n_times = 10,
#'                  offset = NULL,
#'                  b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50)
#' fit_2way$Results
#' #------------------------------------------------------------------------------------
#' #Fit model correctly-specified as the actual simulation model.
#' fit_truemodel <- QRRT(Formula = Ri ~ (x1 + x2) ^ 2 + as.factor(x3),
#'                       Data = Sim_Data, Disperse = 1, beta = NULL, n_times = 1,
#'                    offset = NULL,
#'                    b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50)
#' fit_truemodel$Results
#' #------------------------------------------------------------------------------------
#' # Coefficient of x1 in simulated model is 1; use this fact to illustrate
#' # inclusion of an offset in the model.
#' fit_offset <- QRRT(Formula = Ri ~ x2 + x1:x2 + as.factor(x3),
#'                 Data = Sim_Data, Disperse = 1, beta = NULL, n_times = 2,
#'                 offset = x1,
#'                 b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50)
#' fit_offset$Results
#' #------------------------------------------------------------------------------------
#' round(fit_2way$Results, 3)
#' @references Conteh A, Gavin MC, Solomon J. Quantifying illegal hunting: A novel application of the quantitative randomised response technique.
#' Biological Conservation. 2015;189:16-23.
#'
#' @author Meng Cao, F. Jay Breidt
#' @export QRRT
QRRT = function(Formula, Data, Disperse = 1, beta = NULL, n_times = 1,
                offset = NULL,
                b_distribution = c(6, 7, 4, 2, 2, 1, 1, 1, 1, 25) / 50,
                tolerance = 1e-8){
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

  if(is.null(offset)){
    offset <- rep(0, n)
  }

  #######################
  Ri = model.response(a)
  l1_result = c()

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
    while(abs(l1 - l2) > tolerance){
      l1 <- l2
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
  fr <- rep(NA, length(Ri))
  dfr <- matrix(NA, nrow = length(Ri), ncol = length(beta_1))
  R <- (0:m)
  ###############
  for (i in (1:length(Ri))){
    if (sum(R == Ri[i]) > 0){
      fr[i] <- b_distribution[m + 2] * dpois(Ri[i], lambda[i]) +
        b_distribution[which(R == Ri[i])]
    }else{
      fr[i] = b_distribution[m + 2] * dpois(Ri[i], lambda[i])
    }

    for (j in c(1:length(beta_1))){
      dfr[i, j] <- x[i, j] * 1 / fr[i] * b_distribution[m + 2] *
        (Ri[i] - lambda[i]) * dpois(Ri[i], lambda[i])
    }

  }
  dfr_final <- data.frame(colSums(dfr))
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
                 inverse_fisher_info_report)
  names(result) <- c("Results", "AIC", "Maximized_Log_Likelihood",
                     "Number_Missing_Values",
                    "Sample_Size", "Likelihood_by_Starting_Value", "Score", "Inverse_Fisher_Information")
  return(result)
}
#####End QRRT function.
#--------------------------------------------------------------------------------



