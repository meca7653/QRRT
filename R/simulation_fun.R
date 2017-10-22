#' @name simulation_fun
#' @aliases simulation_fun
#' @title Simulated data
#' @description Simulated data based on design matrix, ball design and true beta
#'
#' @param n sample size
#' @param m max number on green balls
#' @param n_orange number of orange balls
#' @param n_green number of green balls
#' @param x design matrix
#' @param true_beta true parameter of beta
#'
#'
#' @return Data frame contains the observation and design matrix
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
#'
#'
#' @author Meng Cao, F. Jay Breidt
#' @export
simulation_fun = function(n = 250, m = 8, n_orange = 25,
                          n_green = c(6,7,4,2,2,1,1,1,1),
                          x = x, true_beta = true_beta,
                          offset = rep(0, n)){

  balls = c(0:(m+1)) #balls: number on green balls and when balls == m+1 is a orange ball
  n_balls = sum(n_green) + n_orange
  p_balls = c(n_green, n_orange)/n_balls
  lambda<-exp(x%*%true_beta + offset)  #mean of simulated data under true response
  truth<-rpois(n,lambda) #truth response
  Ball = sample(balls, size = n, replace = TRUE, prob = p_balls) #add balls design to true data
  Ri = truth*(Ball == m+1) + Ball*(Ball < m+1) #simulated response Ri: the observation
  result = data.frame(Y = Ri, X = x)
  return(result)

}


