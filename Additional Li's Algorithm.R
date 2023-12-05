# Functions used to compute Li's algorithm.
banach <- function(z_a, z_b, 
                   theta_x, theta_1, 
                   steps = 1000){
  # Compute the minimum of KL-divergence under H_0 for
  # P_0 using Banach fixed-point Theorem. We can
  # compute the minimum using a contraction map T
  # where T(theta_0) = -g(theta_0)*mu + theta_0.
  # This is for the n_a=n_b=p=1 case.
  
  lambda_min <- -theta_x + sum(z_a * theta_1)
  lambda_plus <- theta_x + sum(z_b * theta_1)
  
  D <- 1 + abs(z_a) + abs(z_b)
  mu <- 1 / (D^2 * 2 * (1+1)) # D^2 * 2 * (n_a + n_b)
  
  g <- function(theta_0){
    g1 <- z_a * tanh(sum(z_a * theta_0))
    g2 <- z_b * tanh(sum(z_b * theta_0))
    g3 <- z_a * tanh(-theta_x + sum(z_a * theta_1))
    g4 <- z_b * tanh( theta_x + sum(z_b * theta_1))
    return(g1 + g2 - g3 - g4)
  }
  
  # n_0 = 0
  theta_0 <- 0
  
  for(i in (1:steps)){
    theta_0 <- -(g(theta_0)* mu) + theta_0
  }
  
  return(theta_0)
}

check_domain <- function(z_a, z_b, theta_x, theta_1){
  # Check if the 4-tuple (z_a, z_b, theta_x, theta_1)
  # is in the [5%, 95%]-boundary.
  control <-  abs(-theta_x + z_a * theta_1)
  treatment <- abs(theta_x + z_b * theta_1)
  if(control <= log(sqrt(19)) && treatment <= log(sqrt(19))){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

p <- function(z_a, z_b, theta_x, theta_1, y_a, y_b){
  # Pdf of two group, one sample logistic model.
  p_a_num <- exp((-theta_x + sum(z_a * theta_1)) * (y_a + 1))
  p_a_den <- (exp(2*(-theta_x + sum(z_a * theta_1))) + 1)
  p_a <- p_a_num / p_a_den
  
  p_b_num <- exp(( theta_x + sum(z_b * theta_1)) * (y_b + 1))
  p_b_den <- (exp(2*( theta_x + sum(z_b * theta_1))) + 1)
  p_b <- p_b_num / p_b_den
  return(p_a * p_b)
}

expectation <- function(z_a, z_b, 
                        theta_x, theta_1, 
                        theta_0, theta_min = banach(z_a, z_b, theta_x, theta_1),
                        take_log = TRUE){
  # Compute expectation of P_1/P_0* under P_0.
  gamma <- function(x){
    return(log(cosh(x)))
  }
  
  lambda_min <- -theta_x + sum(z_a * theta_1)
  lambda_plus <- theta_x + sum(z_b * theta_1)
  
  f1 <- gamma(lambda_min  - sum(z_a*theta_min) + sum(z_a*theta_0)) - gamma(sum(z_a * theta_0))
  f2 <- gamma(lambda_plus - sum(z_b*theta_min) + sum(z_b*theta_0)) - gamma(sum(z_b * theta_0))
  a1 <- gamma(sum(z_a * theta_min)) - gamma(lambda_min)
  a2 <- gamma(sum(z_b * theta_min)) - gamma(lambda_plus)
  
  log_ex <- f1 + f2 + a1 + a2
  
  return(if(take_log){log_ex} else{exp(log_ex)})
}