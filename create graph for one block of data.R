# Function that computes theta_0* = theta_min
# using the banach fixed-point theorem.
banach <- function(z_a, z_b, 
                   theta_x, theta_1, 
                   steps = 1000){
  
  lambda_min <- -theta_x + sum(z_a * theta_1)
  lambda_plus <- theta_x + sum(z_b * theta_1)
  
  D <- 1 + abs(z_a) + abs(z_b)
  
  mu <- 1 / (D^2 * 2 * (1+1)) # D^2 * 2 * (n_a + n_b)
  
  g_banach <- function(theta_0){
    g1 <- z_a * tanh(sum(z_a * theta_0))
    g2 <- z_b * tanh(sum(z_b * theta_0))
    g3 <- z_a * tanh(-theta_x + sum(z_a * theta_1))
    g4 <- z_b * tanh( theta_x + sum(z_b * theta_1))
    return(g1 + g2 - g3 - g4)
  }
  
  # n_0 = 0
  theta_0 <- 0
  
  for(i in (1:steps)){
    theta_0 <- -(g_banach(theta_0)* mu) + theta_0
  }
  
  return(theta_0)
}

# Depending on the alternative and covariates
# we can create one block of data.
generate_one_block_outcome <- function(z_a, z_b, theta_x, theta_1){
  # Group a
  p_a <- 1 / (1+exp(-(-theta_x + sum(z_a * theta_1))))
  p_b <- 1 / (1+exp(-( theta_x + sum(z_b * theta_1))))
  y_a <- rbinom(1, 1, p_a)
  y_b <- rbinom(1, 1, p_b)
  return(c(y_a, y_b))
}

# Create the parametrization for logistic
# regression.
r <- function(beta){
  return(exp(beta) / (exp(beta) + exp(-beta)))
}

r_inv <- function(beta){
  return(-1/2 * log(- (beta - 1) / beta))
}

g <- function(mu, z_b, c = 0){
  return(r(c + z_b * r_inv(mu)))
}

# Extra function to compute the probability
P_mu <- function(mu, y_a, gmu, y_b){
  return((mu^y_a) * ((1-mu)^(1-y_a)) * (gmu^y_b) * ((1-gmu)^(1-y_b)))
}

p <- function(x, y){
  return(x^y * (1-x)^(1-y))
}

# Function to keep alpha between zero and one.
cap <- function(alpha, maximum = 1, minimum = 0){
  if(alpha >= maximum){
    return(maximum)
  } 
  if(alpha <= minimum){
    return(minimum)
  }
  return(alpha)
}

# Transform z_b such that it is greater than one.
# We can use this by the symmetry of g.
transform_z_b <- function(z_b){
  if(z_b <= -1){
    z_b <- -z_b
    return(z_b)
  }
  if(z_b >= -1 && z_b <= 0){
    z_b <- -1/z_b
    return(z_b)
  }
  if(z_b >= 0 && z_b <= 1){
    z_b <- 1/z_b
    return(z_b)
  }
  return(z_b)
}

# A seperate function to compute an individual 
# mixture given y_a and y_b.
mixture <- function(alpha, star, corner, y_a, y_b){
  P_a <- (alpha * p(star[1], y_a) + (1 - alpha) * p(corner[1], y_a))
  P_b <- (alpha * p(star[2], y_b) + (1 - alpha) * p(corner[2], y_b))
  return(P_a * P_b)
}


# Once we have computed the mu_star and corner
# probability density we can compute the KL_mixture
# This can also be used for any other mixture.
KL_mix <- function(alpha, alternative, star, corner){
  # alternative, star and corner are two points (mu, g(mu))
  # Define the whole outcome space.
  y_a <- c(0,0,1,1)
  y_b <- c(0,1,0,1)
  
  # To make everything more readable I define the KL
  # here below as one function which we then call
  log_mix <- function(y_a, y_b){
    P_alt <- P_mu(alternative[1], y_a, alternative[2], y_b)
    P_mix <- mixture(alpha, star, corner, y_a, y_b)
    
    return(P_alt * log(P_alt / (P_mix)))
  }
  
  E <- sum(mapply(function(y_a, y_b){log_mix(y_a, y_b)}, 
                  y_a = y_a, y_b = y_b))
  
  return(E)
}



# This function computes the mu_star and determines
# the linear lines. It returns a list. We obtain
# mu_star using Banach's theorem. We now also
# compute the KL-divergence once we found
# mu*, alpha_min and appropriate corner point.
mu_star <- function(z_b, theta_x, theta_1, steps = 10^5){
  D <- 1 + abs(z_b) + 1
  bound <- 1 / (D^2 * 2 * (1+1)) # D^2 * 2 * (n_a + n_b)
  theta_0 <- 0
  
  # We compute mu_star by transforming
  # z_b first, and then transforming
  # mu_star appropriately afterwards.
  # So we take the s-graph where z_b >=1
  # and y_a=y_b=1 as our standard graph.
  trans_z_b <- transform_z_b(z_b)
  
  g_star <- function(theta_0){
    g1 <- trans_z_b * tanh(trans_z_b * theta_0)
    c1 <- trans_z_b
    g2 <- tanh(theta_0)
    c2 <- 1
    return((g1 - c1) + (-g2 +c2))
  }
  
  # Banach Theorem
  for(i in (1:steps)){
    theta_0 <- -(g_star(theta_0) * bound) + theta_0
  }
  
  # Now that we have the maximum angle, we
  # can compute mu_1_star and mu_2_star
  # so the red dots in the graph.
  # We can also get the appropriate corner
  # points.
  mu_1_star <- r(theta_0)
  gmu_1_star <- g(mu = mu_1_star, 
                  z_b = trans_z_b)
  mu_2_star <- 1-mu_1_star
  gmu_2_star <- g(mu = mu_2_star, 
                  z_b = trans_z_b)
  Y_1_corner <- c(0,0)
  Y_2_corner <- c(1,1)
  
  
  # Depending on z_b we need to 
  # transform mu_star's to its appropriate
  # point in [0,1]^2. There are four possible
  # options i.e. graphs.
  # Graph 1 is the graph we take
  # as the general case.
  graph_2 <- (z_b <= -1)
  graph_3 <- (z_b >= 0 && z_b <= 1)
  graph_4 <- (z_b >= -1 && z_b <= 0)
  
  if(graph_2){
    gmu_1_star <- 1 - gmu_1_star
    gmu_2_star <- 1 - gmu_2_star
    Y_1_corner <- c(0,1)
    Y_2_corner <- c(1,0)
  }
  
  if(graph_3){
    mu_1_star <- gmu_1_star
    gmu_1_star <- r(theta_0)
    mu_2_star <- gmu_2_star
    gmu_2_star <- 1 -r(theta_0)
    Y_1_corner <- c(0,0)
    Y_2_corner <- c(1,1)
  }
  
  if(graph_4){
    mu_1_star <- gmu_1_star
    gmu_1_star <- 1-r(theta_0)
    mu_2_star <- gmu_2_star
    gmu_2_star <- r(theta_0)
    Y_1_corner <- c(0,1)
    Y_2_corner <- c(1,0)
  }
  
  # For simplification, we concatenate
  # the coordinates of mu*,mu**
  star_1 <- c(mu_1_star, gmu_1_star)
  star_2 <- c(mu_2_star, gmu_2_star)
  
  # We also do this for the coordinate
  # of the alternative point.
  mu_1 <- r(-theta_x + theta_1)
  gmu_1 <- g(mu = mu_1, 
             z_b = z_b, 
             c = theta_x * (1+z_b))
  
  alternative <- c(mu_1, gmu_1)
  
  # We can now compute the minimum alpha
  # using the abc-formula.
  # Note that it is possible to
  # have P_alt <= P_mu_star so that alpha > 1. 
  #I therefore cap it between 0 and 1.
  
  alpha_min <- function(alternative, star, corner){
    p_1_a <- p(alternative[1], corner[1])
    p_1_b <- p(alternative[2], corner[2])
    p_star_a <- p(star[1], corner[1])
    p_star_b <- p(star[2], corner[2])
    
    a <- 2 * (p_star_a - 1) * (p_star_b - 1)
    b <- 2 * (p_star_a - 1) + 2 * (p_star_b - 1)
    b <- b - p_1_a * (p_star_b - 1) - p_1_b * (p_star_a - 1)
    c <- 2 - p_1_a - p_1_b
    
    plus <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
    minus <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
    
    return(c(cap(plus), cap(minus)))
  }
  
  alpha_1_min <- alpha_min(alternative, 
                           star_1, 
                           Y_1_corner)
  alpha_2_min <- alpha_min(alternative, 
                           star_2, 
                           Y_2_corner)

  alpha_1_pos <- which.min(c(KL_mix(alpha_1_min[1], 
                                    alternative, 
                                    star_1, 
                                    Y_1_corner),
                           KL_mix(alpha_1_min[2], 
                                  alternative, 
                                  star_1, 
                                  Y_1_corner)))
  alpha_2_pos <- which.min(c(KL_mix(alpha_2_min[1], 
                                    alternative, 
                                    star_2, 
                                    Y_2_corner),
                           KL_mix(alpha_2_min[2], 
                                  alternative, 
                                  star_2, 
                                  Y_2_corner)))
  
  
  # Now that we know the best alpha we can
  # compute the KL-divergence of the mixture.
  KL_1_mix <- KL_mix(alpha_1_min[alpha_1_pos],
                     alternative,
                     star_1, 
                     Y_1_corner)
  
  KL_2_mix <- KL_mix(alpha_2_min[alpha_2_pos],
                     alternative,
                     star_2,
                     Y_2_corner)
  
  # The final KL-divergence from the mix.
  KL <- min(KL_1_mix, KL_2_mix)

  # Return the results as a data frame.
  result <- list(KL,
                 KL_1_mix,
                 KL_2_mix,
                 star_1,
                 star_2, 
                 alpha_1_min[alpha_1_pos], 
                 alpha_2_min[alpha_2_pos],
                 Y_1_corner,
                 Y_2_corner)
  names(result) <- c("KL",
                     "KL_1_mix",
                     "KL_2_mix",
                     "star_1",
                     "star_2", 
                     "alpha_1_min", 
                     "alpha_2_min",
                     "Y_1_corner",
                     "Y_2_corner")
  
  return(result)
}

# We are now ready to compute the E-test statistic.
Compute_E_value <- function(z_b, 
                            theta_x, 
                            theta_1, 
                            z_a = 1, steps = 1000){
  # Transform alternative block
  # of data into its parametrization
  mu_1 <- r(-theta_x + theta_1)
  gmu_1 <- g(mu = mu_1, z_b = z_b, 
             c = theta_x * (1+z_b))
  
  # From the given z_b we can compute the convex hull
  # or X-band using the star_mu() function.
  # With the band we can also compute the minimal 
  # alpha in [0,1].
  band <- mu_star(z_b = z_b, 
                  theta_x = theta_x, 
                  theta_1 = theta_1, 
                  steps = 10^5)
  mu_1_star <- band$star_1[1]
  mu_2_star <- band$star_2[1]
  gmu_1_star <- band$star_1[2]
  gmu_2_star <- band$star_2[2]
  alpha_1_min <- band$alpha_1_min
  alpha_2_min <- band$alpha_2_min
  Y_1_corner <- band$Y_1_corner
  Y_2_corner <- band$Y_2_corner

  # We are now ready to compute the KL mixtures.
  alternative <- c(mu_1, gmu_1)
  star_1 <- c(mu_1_star, gmu_1_star)
  star_2 <- c(mu_2_star, gmu_2_star)
  
  # I am going to select a few alphas and choose estimate
  # the minimum alpha we need.
  alpha <- (1:100) / 100
  KL_alpha_1 <- mapply(function(alpha){KL_mix(alpha, alternative, 
                                            star_1, Y_1_corner)}, 
                     alpha)
  alpha_1_min_est <- alpha[which.min(KL_alpha_1)]
  
  KL_alpha_2 <- mapply(function(alpha){KL_mix(alpha, alternative, 
                                              star_2, Y_2_corner)}, 
                       alpha)
  alpha_2_min_est <- alpha[which.min(KL_alpha_2)]
  
  KL_1_est <- KL_mix(alpha_1_min_est, alternative,
                     star_1, Y_1_corner)
  KL_2_est <- KL_mix(alpha_2_min_est, alternative,
                     star_2, Y_2_corner)
  
  # This chooses one of the two alpha's it computed
  # from the quadratic formula.
  KL_mu_1_star <- band$KL_1_mix
  KL_mu_2_star <- band$KL_2_mix

  # We also need to compute the KL-divergence for theta_0* = theta_min.
  # For this alpha = 1 and mu_0 = r(theta_0*).
  theta_min <- banach(1, z_b, 
                      theta_x, theta_1, 
                      steps = 10^5)
  mu_min <- r(theta_min)
  gmu_min <- g(mu_min, z_b)
  KL_min <- KL_mix(1, alternative, 
                   c(mu_min, gmu_min), c(1,1))
  
  print(alpha_1_min)
  print(KL_alpha_1[which(KL_alpha_1 == alpha_1_min)])
  print(alpha_2_min)
  print(KL_alpha_2[which(KL_alpha_2 == alpha_2_min)])
  
  
  # Now we have computed the KL divergence and are able to choose
  # the one that is the lowest. With this we can compute
  # each LR such that we have an e-value. For this we need
  # the outcome y = (y_a, y_b). Generate outcome y:
  y <- generate_one_block_outcome(z_a, z_b, theta_x, theta_1)
  #y <- c(1,1)

  # Lastly, we compute the individual LR statistics.
  LR_star_1 <- P_1 / mixture(alpha_1_min, star_1, 
                             Y_1_corner, y[1], y[2])
  LR_star_2 <- P_1 / mixture(alpha_2_min, star_2, 
                             Y_2_corner, y[1], y[2])
  LR_min <- P_1 / P_mu_min
  LR_est_1 <- P_1 / mixture(alpha_1_min_est, star_1, 
                            Y_1_corner, y[1], y[2])
  LR_est_2 <- P_1 / mixture(alpha_2_min_est, star_2, 
                            Y_2_corner, y[1], y[2]) 
  

  # The minium KL-divergence is chosen as the one that gives
  # the LR.
  minimum_KL <- which.min(c(KL_mu_1_star,KL_mu_2_star, 
                            KL_min, KL_1_est, KL_2_est))
  
  # Printing jobs to see that we do it correct.
  print("--------KL mixtures-----------")
  print(c("KL_mu_1_star","KL_mu_2_star",
          "KL_min","KL_1_est", "KL_2_est"))
  print(c(KL_mu_1_star,KL_mu_2_star, 
          KL_min,  KL_1_est, KL_2_est))
  print(c("KL_mu_1_star","KL_mu_2_star",
          "KL_min","KL_1_est", "KL_2_est")[minimum_KL])
  print("-------------LR-------------")
  print(c(LR_star_1,LR_star_2, 
          LR_min, LR_est_1, LR_est_2))

  E_value <- c(LR_star_1,LR_star_2, LR_min, 
               LR_est_1, LR_est_2)[minimum_KL]
  
  #------------------ PLOTTING -----------------------------
  # To view everything I plot g with its appropriate bands.
  mu_0 <- r((-10000:10000) / 100)
  gmu_0 <- g(mu = mu_0, z_b = z_b)
  t <- (0:500) / 500
  
  plot(mu_0^y[1] * (1-mu_0)^(1-y[1]), 
       gmu_0^y[2] * (1-gmu_0)^(1-y[2]), 
       type = 'l', 
       col = 'black', 
       xlab = 'mu = r(beta)', 
       ylab = 'g(mu)', 
       main = paste('z_b =', round(z_b,2), 
                    ', y_a = ',y[1], 
                    ', y_b =',y[2],
                    '\n theta_x =', round(theta_x,2),
                    ', theta_1 =', round(theta_1,2),
                    '\n E_value =',round(E_value,2)))
  
  # In red, the boundary that determines
  # the convex hull.
  points(p(mu_1_star, y[1]), 
         p(gmu_1_star, y[2]), 
         col = 'red', pch = 19)
  lines(t*p(mu_1_star, 
            y[1]) + (1-t)*p(Y_1_corner[1], 
                            y[1]), 
        t*p(gmu_1_star, 
            y[2]) + (1-t)*p(Y_1_corner[2], 
                            y[2]), 
        col = 'red')
  points(p(mu_2_star, y[1]), 
         p(gmu_2_star, y[2]), 
         col = 'red', pch = 19)
  lines(t*p(mu_2_star, 
            y[1]) + (1-t)*p(Y_2_corner[1], 
                                       y[1]), 
        t*p(gmu_2_star, 
            y[2]) + (1-t)*p(Y_2_corner[2], 
                                        y[2]), 
        col = 'red')
  
  # In purpble, the minimum obtained
  # from minimizing the KL-divergence
  # for H_0.
  points(p(mu_min, y[1]), 
         p(gmu_min, y[2]), 
         col = 'purple', 
         pch = 19)
  
  # In orange, the estimated alpha's.
  points(p(alpha_1_min_est * mu_1_star + (1-alpha_1_min_est) * Y_1_corner[1], 
           y[1]), 
         p(alpha_1_min_est * gmu_1_star + (1-alpha_1_min_est) * Y_1_corner[2], 
           y[2]), 
         col = 'orange', pch = 19)
  points(p(alpha_2_min_est * mu_2_star + (1-alpha_2_min_est) * Y_2_corner[1], 
           y[1]), 
         p(alpha_2_min_est * gmu_2_star + (1-alpha_2_min_est) * Y_2_corner[2], 
           y[2]), 
         col = 'orange', 
         pch = 18, cex = 1.5)
  
  # Above the estimated alpha's, I plot 
  # the correct alpha's in green.
  points(p(alpha_1_min * mu_1_star + (1-alpha_1_min) * Y_1_corner[1], 
           y[1]), 
         p(alpha_1_min * gmu_1_star + (1-alpha_1_min) * Y_1_corner[2], 
           y[2]), 
         col = 'darkgreen', 
         pch = 19)
  points(p(alpha_2_min * mu_2_star + (1-alpha_2_min) * Y_2_corner[1], 
           y[1]), 
         p(alpha_2_min * gmu_2_star + (1-alpha_2_min) * Y_2_corner[2], 
           y[2]), 
         col = 'darkgreen', 
         pch = 18, 
         cex = 1.5)
  
  # The alternative
  points(p(mu_1, y[1]), 
         p(gmu_1, y[2]), 
         col = 'blue')
  
  return(E_value)
}


# Pick some random z_b, theta_x and theta_1
z_b <- rnorm(1, sd = 4)
theta_x <- rnorm(1)
theta_1 <- rnorm(1)

# Compute the E-test statistic.
Compute_E_value(z_b, theta_x, theta_1)

#Compute_E_value(z_b = 2.5, theta_x = .1, theta_1 = .1)

