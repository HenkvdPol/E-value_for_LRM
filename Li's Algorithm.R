source("~/Desktop/Additional Li's Algorithm.R")
set.seed(123)

Li <- function(z_a, z_b, theta_x, theta_1, 
               max_steps = 10, 
               size_beta = 100, 
               size_alpha = 100, 
               threshold = 0){
  # Compute the mixture minimum using Li's algorithm.
  
  # First check if we have a 4-tuple that is within our boundary setting.
  if(check_domain(z_a, z_b, theta_x, theta_1) == FALSE){
    return(print("Out of bound."))
  }
  
  # Parameter settings.
  theta_min <- banach(z_a, z_b, theta_x, theta_1) # Initial minimum of the KL-divergence.
  y_a = c(-1,1,-1,1) # Input for the function.
  y_b = c(-1,1,1,-1) # So that it stays consistent throughout the function.
  beta <- rnorm(5*10^5, theta_min, 1/10) # Initial beta's to check.
  best_beta <- rep(0, max_steps) # parameter to save the minimum beta.
  best_beta[1] <- theta_min # Make P_1 = P_theta_min.
  alpha <- (1:size_alpha) / size_alpha # Hundred points in [0,1].
  best_alpha <- rep(0, max_steps) # parameter to save the minimum alpha.
  
  # about 20 twenty points of theta_0 we check if it is an E-value.
  theta_0 <- c((-5:5)/5, -5,5,100,-100, theta_min + 10^-4, theta_min - 10^-4)
  # Only need to compute this once.
  P_theta_1 <- mapply(function(y_a, y_b){p(z_a, z_b, theta_x, theta_1, y_a, y_b)}, 
                      y_a = y_a, y_b = y_b)
  # Keep track of the best iteration.
  best = c(0,100)
  # Look how max(E_theta_0) has decreased in each iteration.
  E_max <- rep(0, max_steps)
  
  # Recursion formula to compute the mixtures.
  P_m <- function(y_a, y_b, best_alpha, best_beta, i){
    if(i == 0){
      return(0)
    }
    else{
      i = i - 1
      # Some print statements that are sometimes useful.
      #print(paste("alpha:", best_alpha[i+1]))
      #print(paste("beta:", best_beta[i+1]))
      
      a <- best_alpha[i+1]
      P_min <- P_m(y_a, y_b, best_alpha, best_beta, i)
      P_b <- p(z_a, z_b, 0, best_beta[i+1], y_a, y_b)
      # Mixture
      return(a * P_min + (1 - a) * P_b)
    }
  }
  
  # D(P_theta_1 || alpha * P_(m-1) + (1-alpha) * P_beta)
  mixed_KL <- function(alpha, beta, Pm){
    # I can leave out p_(theta_x, theta_1) in the log and find the minimum.
    mix <- function(y_a, y_b, j){
      P_b <- p(z_a, z_b, 0, beta, y_a, y_b)
      return(- P_theta_1[j] * log((alpha * Pm[j] + (1-alpha) * P_b)))
    }
    # We now have n_a = n_b = 1.
    result <- mapply(mix, y_a = y_a, y_b = y_b, j = c(1,2,3,4))
    return(rowSums(result))
  }
  
  # Compute the expectation: E_theta_0[P_theta_1 / P_(m)]
  E_statistic <- function(theta_0, Pm){
    ex <- function(y_a, y_b, j){
      return(p(z_a, z_b, 0, theta_0, y_a, y_b) * P_theta_1[j] / Pm[j])
    }
    E_result <- mapply(ex, y_a = y_a, y_b = y_b, j = c(1,2,3,4))
    return(log(sum(E_result)))
  }  
  
  for(i in 1:max_steps){
    #if(i %% i == 0){print(paste("step:",i))}
    # Before we continue to the next iteration, we first check if the current is correct.
    Pm <- mapply(function(y_a, y_b){P_m(y_a, y_b, best_alpha, best_beta, i)}, 
                 y_a = y_a, y_b = y_b) # only compute this once.
    E_theta_0 <- mapply(function(theta_0){E_statistic(theta_0, Pm = Pm)}, 
                        theta_0 = theta_0)
    
    # Check if we have improved in this iteration compared to the previous ones.
    if(max(E_theta_0) <= best[2]){
      best[1] <- i
      best[2] <- max(E_theta_0)
    }
    
    # check if we are done.
    if(all(E_theta_0 <= threshold)){
      print("ladys and gentlemen...")
      break
    }
    
    # Quit if we have reached max_step to avoid too many computations.
    if(i / max_steps == 1){next}
    
    #print(paste("max E_theta_0:", max(E_theta_0)))
    E_max[i] <- max(E_theta_0)
    #print("--------------------")
    # Step one: compute the minimum KL-divergence mixture
    # by taking the minimum of the KL mixture in a grid.
    grid <- sapply(beta, function(beta){mixed_KL(alpha = alpha, beta, Pm)})
    minimum <- which.min(grid)
    
    # Step two: Retrieve the optimal alpha and beta.
    best_alpha[i+1] <- alpha[as.integer(minimum %% length(alpha)) + 1]
    best_beta[i+1] <- beta[as.integer(floor(minimum / length(alpha))) + 1]
    
    #print(paste("best alpha in step", i, "is", best_alpha[i+1]))
    #print(paste("best beta in step", i, "is", best_beta[i+1]))
    
    # Step three: Update new beta. Because of the iteration step, 
    # we can not compute a new theta_min. I therefore pick randomly
    # new points around theta_min where I increase the sd of
    # the new points. So that if we can not find the better beta 
    # in the first steps. We can pick a new one that is further away.
    beta <- rnorm(size_beta, mean = theta_min, sd = i*sqrt(i))
  }
  print("Results")
  print("Best alpha choices.")
  print(best_alpha)
  print("Best Beta choices.")
  print(best_beta)
  print(paste("best iteration:", best))
  print(paste("original", max(mapply(function(theta_0){expectation(z_a, z_b, 
                                                                   theta_x, theta_1, 
                                                                   theta_0)}, 
                                     theta_0 = (-100:100) / 10))))
  
  # At the end we are interested in how our expectation now looks like, compared to 
  # the original one.
  Pm <- mapply(function(y_a, y_b){P_m(y_a, y_b, 
                                      best_alpha, best_beta, 
                                      best[1])}, 
               y_a = y_a, y_b = y_b)
  plot((-100:100)/10, mapply(function(theta_0){E_statistic(theta_0, Pm = Pm)}, 
                             theta_0 = (-100:100) / 10), 
       ylab = "log(Expectation)", xlab = "Theta_0", type = 'l', col = 'darkblue')
  
  # The original minimum.
  lines((-100:100)/10, mapply(function(theta_0){expectation(z_a, z_b, 
                                                            theta_x, theta_1, 
                                                            theta_0, theta_min)}, 
                              theta_0 = (-100:100)/10), 
        col = 'orange')
  
  # See how many times we choose alpha's or beta's to get an indication of where better e-values are.
  hist(best_alpha, breaks = 50)
  hist(best_beta, breaks = 50)
  # Look how the maximum of the expectation is in each iteration.
  plot((1:max_steps), E_max, xlab = 'iteration', ylab = 'max Expectation', type = 'l')
  # Did Li's algorithm succeed?
  return(all(mapply(function(theta_0){E_statistic(theta_0, Pm = Pm)}, 
                    theta_0 = (-100:100) / 10) <= threshold))
}

Li(1,2,1/2, -1/2, max_steps = 10, size_beta = 10, size_alpha = 100)