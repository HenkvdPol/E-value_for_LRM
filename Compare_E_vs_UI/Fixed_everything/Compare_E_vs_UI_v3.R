# Some functions we use throughout to make everything
# more readable. (like finding theta_0*)
source("additional_E_vs_UI_functions.R")

# Set seed to get comparative results
set.seed(123)

# This is the same code as in V2. However, we now
# make the loop a function such that, when we see
# a statistic that is greater then 20. We stop 
# looping through it. We then check for which values
# we stop.

Sample_data <- function(stopping_value = 20, maximum_loops = 100){
  # Create data arrays including two dummy arrays
  E_results_2 <- rep(NA, maximum_loops)
  E_i_results <- rep(NA, maximum_loops)
  UI_results_2 <- rep(NA, maximum_loops)
  UI_j_results <- rep(NA, maximum_loops)
  
  # We start with initializing some data aalready.
  # Such that we do not get 0,1-probabilities.
  Y <- c(c(0,0), c(0,1), c(1,0), c(1,1))
  begin <- 4
  z_a <- rep(1,4) # rnorm(4,1,2) # 
  z_b <- rep(2,4) # rnorm(4,1,2) # 
  true_theta_x <- rep(1/2, 4) # rnorm(4,1,1/4) # 
  true_theta_1 <- rep(1/2, 4) # rnorm(4,1,1/4) # 
  

  for(i in 1:maximum_loops){
    #print(paste("Step", i))
    # For each iteration, we need to do a couple of things:
    # (i) We need to estimate the MLE of the past seen data under
    # the alternative.
    alt_old <- estimate_coefficients(Y, # (i - 1)-data points
                                     z_a, # + begin
                                     z_b)
    #alt_old <- c(1/2, 1/2)
    #print(paste("est theta_x:", alt_old[1]))
    #print(paste("est theta_1:", alt_old[2]))
  
    # (ii) collect new data.
    z_a_new <- 1 # rnorm(1,1,2) # 
    z_b_new <- 2 # rnorm(1,1,2) # 
    true_theta_x_new = 1/2 # rnorm(1,1,1/4) # 
    true_theta_1_new = 1/2 # rnorm(1,1,1/4) # 
    Y_new <- generate_one_block_outcome(z_a_new, 
                                        z_b_new, 
                                        true_theta_x_new, 
                                        true_theta_1_new)
    
    
    # (iii) Update the data parameters.
    #true_theta_x <- c(true_theta_x, true_theta_x_new)
    #true_theta_1 <- c(true_theta_1, true_theta_1_new)
    z_a <- c(z_a, z_a_new)
    z_b <- c(z_b, z_b_new)
    Y <- c(Y, Y_new)
    #est_alternatives <- c(est_alternatives, alt_old)
    
    
    # (iv) estimate the MLE of all the seen data under the null.
    theta_null <- estimate_coefficients(Y[(2*4 + 1):length(Y)], # i data points
                                        z_a[(4 + 1):length(z_a)], # without begin.
                                        z_b[(4 + 1):length(z_b)], 
                                        null = 0)
    #print(paste("Theta_null :", theta_null))
    
    # (v) Compute E-value and UI-value.
    E_i_results[i] <- Compute_E_i(z_a_new, 
                                  z_b_new, 
                                  alt_old[1], 
                                  alt_old[2], 
                                  Y_new)
    E_results_2[i] <- prod(E_i_results, na.rm = TRUE)
    #print(E_i_results)
    
    
    # We need to loop for all the probabilities in each step
    # for UI as alt_old and theta_null changes when we retrieve
    # more data. 
    for(j in begin:(length(z_b)-1)){ # j = 1,...,i
      MLE_mu_null <- parametrized_null(z_a[j + 1], 
                                       z_b[j + 1], 
                                       theta_null)
      
      P_null <- P_mu(MLE_mu_null[1], Y[j*2 + 1],
                     MLE_mu_null[2], Y[j*2 + 2])
      
      MLE_mu_alt <- parametrized_alternative(z_a[j + 1], 
                                             z_b[j + 1], 
                                             alt_old[1], 
                                             alt_old[2])
      P_1 <- P_mu(MLE_mu_alt[1], Y[j*2 + 1],
                  MLE_mu_alt[2], Y[j*2 + 2])
      
      UI_j_results[j] <- P_1 / P_null
    }
    #print(UI_j_results)
    UI_results_2[i] <- prod(UI_j_results, na.rm = TRUE)
    
    # If one of the results is greater than the stopping value.
    # We stop the loop and return.
    if(E_results_2[i] >= stopping_value){
      print(paste('E_i :', E_results_2[i], 'at step', i))
      return(c('E',i))
    }
    if(UI_results_2[i] >= stopping_value){
      print(paste('UI_i:', UI_results_2[i], 'at step', i))
      return(c('UI', i))
    }
    
    if(all(c(E_results_2[i], UI_results_2[i]) >= stopping_value)){
      print(paste('at step', i))
      print(paste('E_i :', E_results_2[i]))
      print(paste('UI_i:', UI_results_2[i]))
      return(c('UI and E', i))
    }
    
    # After each iteration, we reset the UI dummy array.
    UI_j_results <- rep(NA, maximum_loops)
  }
  # If after everything, the loop does not hold any significant E-value. We quit.
  return(c("No significance", i))
}

start.time <- Sys.time()
Hello <- replicate(1000, Sample_data(stopping_value = 20, maximum_loops = 100))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
write.csv(Hello)


