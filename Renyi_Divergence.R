check_domain <- function(z_a, z_b, theta_x, theta_1){
  # Control group and Treatment Group
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

p_i <- function(y_i, z_i, x, theta_x, theta_1){
  p_num <- exp((x * theta_x + sum(z_i * theta_1)) * (y_i + 1))
  p_den <- (exp(2*(x * theta_x + sum(z_i * theta_1))) + 1)
  return(p_num / p_den)
}

Renyi_Banach <- function(z_a, z_b, 
                         theta_x, theta_1, 
                         alpha, beta, 
                         steps = 1000){
  
  lambda_min <- -theta_x + sum(z_a * theta_1)
  lambda_plus <- theta_x + sum(z_b * theta_1)
  
  D <- 1 + abs(z_a) + abs(z_b) + (1+abs(alpha)) + (1+abs(beta))
  
  mu <- 1 / (D^2 * 2 * (1+1)) # D^2 * 2 * (n_a + n_b)
  
  g <- function(theta_0){
    g1 <- z_a * tanh(sum(z_a * theta_0))
    g2 <- z_b * tanh(sum(z_b * theta_0))
    g3 <- z_a * tanh(lambda_min * alpha + z_a * theta_0 * (1- alpha))
    g4 <- z_b * tanh(lambda_plus * beta + z_b * theta_0 * (1- beta))
    return(g1 + g2 - g3 - g4)
  }
  
  # n_0 = 0
  theta_0 <- 0
  
  for(i in (1:steps)){
    theta_0 <- -(g(theta_0)* mu) + theta_0
  }
  
  return(theta_0)
}

M_value <- function(y_a, y_b, 
                    z_a, z_b, 
                    theta_x, theta_1, 
                    alpha, beta, 
                    nsteps = 1000,
                    theta_min = Renyi_Banach(z_a, 
                                             z_b, 
                                             theta_x, 
                                             theta_1, 
                                             alpha, 
                                             beta, 
                                             nsteps)){
  y <- c(-1,1)
  
  M_A_num <- (p_i(y_a, z_a, 
                  -1, theta_x, 
                  theta_1) / p_i(y_a, z_a, 
                                 -1, 0, 
                                 theta_min))^alpha
  M_A_den <- sum(mapply(function(y){p_i(y, 
                                        z_a, 
                                        -1, 
                                        theta_x, 
                                        theta_1)^alpha * p_i(y, 
                                                             z_a, 
                                                             -1, 
                                                             0, 
                                                             theta_min)^(1-alpha)}, 
                        y = y))
  
  M_B_num <- (p_i(y_b, z_b, 
                  1, theta_x, 
                  theta_1) / p_i(y_b, z_b, 
                                 1, 0, 
                                 theta_min))^beta
  M_B_den <- sum(mapply(function(y){p_i(y, 
                                        z_b, 
                                        1, 
                                        theta_x, 
                                        theta_1)^beta * p_i(y, 
                                                            z_b, 
                                                            1, 
                                                            0, 
                                                            theta_min)^(1-beta)}, 
                        y = y))
  
  return((M_A_num / M_A_den) * (M_B_num / M_B_den))
}

Expected_M_value <- function(theta_0, 
                             z_a, z_b, 
                             theta_x, theta_1, 
                             alpha, beta, 
                             nsteps = 1000,
                             theta_min = Renyi_Banach(z_a, 
                                                      z_b, 
                                                      theta_x, 
                                                      theta_1, 
                                                      alpha, 
                                                      beta, 
                                                      nsteps)){
  y_a <- c(-1,-1,1,1)
  y_b <- c(-1, 1,-1,1)
  
  Expectation <- mapply(function(y_a, y_b){p(z_a, 
                                             z_b, 
                                             0, 
                                             theta_0, 
                                             y_a, 
                                             y_b) * M_value(y_a,
                                                            y_b, 
                                                            z_a,
                                                            z_b, 
                                                            theta_x, 
                                                            theta_1, 
                                                            alpha, 
                                                            beta, 
                                                            theta_min = theta_min)}, 
                        y_a = y_a, y_b = y_b)
  return(sum(Expectation))
}

Renyi_Stoplight <- function(z_a, z_b, alpha, beta, nsteps = 1000, threshold = 1){
  results <- matrix(nrow = length((-75:75) / 50), ncol = length((-75:75) / 50))
  minimum_results <- results
  
  for(x in 1:151){
    for(y in 1:101){
      theta_x <- (x - 75) / 50
      theta_1 <- (y - 50) / 50
      if(check_domain(z_a, z_b, theta_x, theta_1) == FALSE){next}
      
      theta_min <- Renyi_Banach(z_a, z_b, 
                                theta_x, theta_1, 
                                alpha, beta, 
                                steps = nsteps)
      
      minimum_results[x,y] <- theta_min
      
      t_0 <- c((-5:5)/5, -5,5,100,-100, theta_min + 10^-4, theta_min - 10^-4)
      ex <- mapply(function(t_0){Expected_M_value(t_0, 
                                                  z_a, 
                                                  z_b, 
                                                  theta_x, 
                                                  theta_1, 
                                                  alpha, 
                                                  beta, 
                                                  theta_min = theta_min)}, 
                   t_0 = t_0)
      
      results[x,y] <- max(ex)
    }
  }
  print("Minimum values")
  write.csv(minimum_results)
  print('--------------------------------------- RESULTS -----------------------------------------------------')
  print('--------------------------------------- RESULTS -----------------------------------------------------')
  print('--------------------------------------- RESULTS -----------------------------------------------------')
  write.csv(results)
}

start <- proc.time()
Renyi_Stoplight(1,2,1/4,1/4, nsteps = 10^5)
print(proc.time() - start)


