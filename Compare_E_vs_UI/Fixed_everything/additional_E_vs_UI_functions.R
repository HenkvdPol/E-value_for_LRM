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

check_domain <- function(z_a, z_b, theta_x, theta_1){
  # Control group and Treatment Group
  control <-  abs(-theta_x + z_a * theta_1)
  treatment <- abs(theta_x + z_b * theta_1)
  
  # Linear condition to hold for [5%, 95%]-interval.
  if(control <= log(sqrt(19)) && treatment <= log(sqrt(19))){
    return(TRUE)
  }
  else{
    return(FALSE)
  }
}

# Create the parametrization
r <- function(beta){
  result <- exp(beta) / (exp(beta) + exp(-beta))
  # Check if we do not get NaN values. We can only
  # get NaN values if we have inf/inf-situation.
  # In this situation r(beta) = 1.
  if(is.na(result)){return(1)}
  # Otherwise, just return the computed value.
  return(result)
}

r_inv <- function(beta){
  return(-1/2 * log(- (beta - 1) / beta))
}

g <- function(mu, z_b, c = 0){
  return(r(c + z_b * r_inv(mu)))
}

parametrized_null <- function(z_a, z_b, 
                              theta_0){
  mu_0 <- r(z_a * theta_0)
  # This is to prevent dividing by zero.
  if(mu_0 == 1){mu_0 <- .99999999999}
  gmu_0 <- r((z_b / z_a) * r_inv(mu_0))
  return(c(mu_0, gmu_0))
}

parametrized_alternative <- function(z_a, z_b, 
                                     theta_x, 
                                     theta_1){
  mu_1 <- r(-theta_x + z_a * theta_1)
  # This is to prevent dividing by zero.
  if(mu_1 == 1){mu_1 <- .99999999999}
  
  z <-  z_b / z_a
  gmu_1 <- r(theta_x * (1 + z) + z * r_inv(mu_1))
  return(c(mu_1, gmu_1))
}

# Extra function to compute the probability
P_mu <- function(mu, y_a, gmu, y_b){
  p_1 <- mu^y_a 
  p_2 <- (1-mu)^(1-y_a)
  p_3 <- gmu^y_b
  p_4 <- (1-gmu)^(1-y_b)
  return(p_1 * p_2 * p_3  * p_4)
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

# A seperate function to compute an individual 
# mixture given y_a and y_b.
mixture <- function(alpha, star, corner, y_a, y_b){
  # star and corner are two-dimensional arrays.
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
    P_alt <- P_mu(alternative[1], y_a, 
                  alternative[2], y_b)
    P_mix <- mixture(alpha, star, 
                     corner, y_a, y_b)
    
    # Extra condition so that we do not produce
    # NaN when doing mixtures.
    if(P_alt == 0){return(0)}
    
    return(P_alt * log(P_alt / (P_mix)))
  }
  
  E <- sum(mapply(function(y_a, y_b){log_mix(y_a, y_b)}, 
                  y_a = y_a, y_b = y_b))
  
  return(E)
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

# Depending on the alternative and covariates
# we can create one block of data. I do this
# here without the parametrization.
generate_one_block_outcome <- function(z_a, z_b, theta_x, theta_1){
  # Group a
  p_a <- 1 / (1+exp(-2 * (-theta_x + z_a * theta_1)))
  p_b <- 1 / (1+exp(-2 * ( theta_x + z_b * theta_1)))
  y_a <- rbinom(1, 1, p_a)
  y_b <- rbinom(1, 1, p_b)
  return(c(y_a, y_b))
}

# This function computes the P_mix function by
# first computing mu_star that determines
# the linear lines of conv(H_0). It returns a list.
# mu_star is obtained using Banach's theorem.
P_mixture <- function(z_a, z_b, 
                      theta_x, 
                      theta_1, 
                      y_a, 
                      y_b,
                      nsteps = 10^5){
  # For one data block, we can transform
  # z_b by z_a such that z_a = 1 and z_b = z_b/z_a.
  z <- z_b / z_a
  
  # Parameters for the contraction map.
  D <- 1 + abs(z_b) + 1
  bound <- 1 / (D^2 * 2 * (1+1)) # D^2 * 2 * (n_a + n_b)
  theta_0 <- 0
  
  # We compute mu_star by transforming
  # z_b first, and then transforming
  # mu_star appropriately afterwards.
  # So we take the s-graph where z_b/z_a >=1
  # and y_a=y_b=1 as our standard graph.
  trans_z_b <- transform_z_b(z)
  
  g_star <- function(theta_0){
    g1 <- trans_z_b * tanh(trans_z_b * theta_0)
    c1 <- trans_z_b
    g2 <- tanh(theta_0)
    c2 <- 1
    return((g1 - c1) + (-g2 +c2))
  }
  
  # Banach Theorem
  for(i in (1:nsteps)){
    theta_0 <- -(g_star(theta_0) * bound) + theta_0
  }
  
  # Now that we have the maximum angle, we
  # can compute mu_1_star and mu_2_star
  # so the red dots in the graph.
  # We can also get the appropriate corner
  # points.
  star_1 <- parametrized_null(1, trans_z_b, theta_0)
  mu_1_star <- r(theta_0)
  gmu_1_star <- g(mu = mu_1_star, 
                  z_b = trans_z_b)
  star_2 <- parametrized_null(1, trans_z_b, theta_0)
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
  graph_2 <- (z <= -1)
  graph_3 <- (z >=  0 && z <= 1)
  graph_4 <- (z >= -1 && z <= 0)
  
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
  # of the alternative point. Here we do need
  # the distinct z_a, z_b. Since if z_a \neq 1.
  # Then the alternative point changes with it.
  alternative <- parametrized_alternative(z_a, 
                                          z_b, 
                                          theta_x, 
                                          theta_1)
  # We can now compute the minimum alpha
  # using the abc-formula.
  # Note that it is possible to
  # have P_alt <= P_mu_star so that alpha > 1. 
  # I therefore cap it between 0 and 1.
  # Note that the discriminant is always non-negative.
  alpha_min <- function(alternative, star, corner){
    p_1_a <- p(alternative[1], corner[1])
    p_1_b <- p(alternative[2], corner[2])
    p_star_a <- p(star[1], corner[1])
    p_star_b <- p(star[2], corner[2])
    
    a <- 2 * (p_star_a - 1) * (p_star_b - 1)
    b <- (2 - p_1_b) * (p_star_a - 1) + (2 - p_1_a) * (p_star_b - 1)
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
  
  # Check which alpha returns the best KL-mixture.
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
  KL_mix_1 <- KL_mix(alpha_1_min[alpha_1_pos],
                     alternative,
                     star_1, 
                     Y_1_corner)
  
  KL_mix_2 <- KL_mix(alpha_2_min[alpha_2_pos],
                     alternative,
                     star_2,
                     Y_2_corner)
  
  # We want to compare it to the minimum
  # of the KL-divergence under the null.
  # Here we also need the given z_a, z_b.
  theta_min <- banach(z_a, 
                      z_b, 
                      theta_x, 
                      theta_1, 
                      steps = nsteps)
  
  mu_null <- parametrized_null(z_a, 
                               z_b, 
                               theta_min)
  
  KL_min <- KL_mix(1, 
                   alternative, 
                   mu_null, 
                   c(1,1))
  
  #  One of the three mixtures has the lowest KL-divergence.
  #  We want to use the lowest mixture probability.
  P_mix <- c(mixture(alpha_1_min[alpha_1_pos], # First mixture
                     star_1, 
                     Y_1_corner, 
                     y_a, y_b),
             mixture(alpha_2_min[alpha_2_pos], # Second mixture
                     star_2, 
                     Y_2_corner, 
                     y_a, y_b),
             mixture(1, 
                     mu_null, # minimum KL-divergence
                     c(1,1), 
                     y_a, y_b))[which.min(c(KL_mix_1,
                                            KL_mix_2,
                                            KL_min))]
  
  # Return the results as a data frame.
  result <- list(P_mix,
                 star_1,
                 star_2, 
                 alpha_1_min[alpha_1_pos], 
                 alpha_2_min[alpha_2_pos],
                 Y_1_corner,
                 Y_2_corner)
  names(result) <- c("P_mix",
                     "star_1",
                     "star_2", 
                     "alpha_1_min", 
                     "alpha_2_min",
                     "Y_1_corner",
                     "Y_2_corner")
  
  return(result)
}

# We also require to compute if the alternative is inside
# or outside the band that is determined by mu_star and g.
inside_band <- function(z_a, 
                        z_b, 
                        theta_x, 
                        theta_1,
                        band = P_mixture(z_a, 
                                         z_b, 
                                         theta_x, 
                                         theta_1, 
                                         y_a,
                                         y_b)){
  # Alternative point in [0,1]^2
  mu_1 <- r(-theta_x + z_a * theta_1)
  # This is to prevent dividing by zero.
  if(mu_1 == 1){mu_1 <- .99999999999}
  # We also adjust it to the correct outcome
  z <-  z_b / z_a
  gmu_1 <- r(theta_x * (1 + z) + z * r_inv(mu_1))
  
  # We have two points such that we can determine the linear line
  # between these points.
  left_band <- function(mu){
    # Determine slope and intercept for the linear line.
    slope <- (band$star_1[2] - band$Y_1_corner[2]) / (band$star_1[1] - band$Y_1_corner[1])
    intercept <- band$Y_1_corner[2] - band$Y_1_corner[1] * slope
    
    # The linear line is only defined for [mu_star, 1] or [0,mu_star]
    upper <- max(c(band$star_1[1], band$Y_1_corner[1]))
    lower <- min(c(band$star_1[1], band$Y_1_corner[1]))
    if(mu >=  upper || mu <= lower){
      # If it is outside, we only look at g(mu).
      return(r((z_b / z_a) * r_inv(mu)))
    }
    else{
      return(intercept + slope * mu)
    }}
  
  # Now we do the same computation for the other side of the band.
  right_band <- function(mu){
    slope <- (band$star_2[2] - band$Y_2_corner[2]) / (band$star_2[1] - band$Y_2_corner[1])
    intercept <- band$Y_2_corner[2] - band$Y_2_corner[1] * slope
    
    lower <- min(c(band$star_2[1], band$Y_2_corner[1]))
    upper <- max(c(band$star_2[1], band$Y_2_corner[1]))
    if(mu <= lower || mu >= upper){
      return(r((z_b / z_a) * r_inv(mu)))
    }
    else{
      return(intercept + slope * mu)
    }}
  
  # The alternative is now inside the band if at mu, it is in between
  # one of the three functions. We also here need to adjust for the
  # outcome block Y=(y_a,y_b). We don't need Y_a necessarily.
  g_left <- left_band(mu_1)
  g_right <- right_band(mu_1)
  g_middle <- r((z_b / z_a) * r_inv(mu_1))
  
  if(gmu_1 >= min(c(g_left, 
                    g_right, 
                    g_middle)) && gmu_1 <= max(c(g_left, 
                                                 g_right, 
                                                 g_middle))){
    return(TRUE)
  }
  # Otherwise
  return(FALSE)
}


estimate_coefficients <- function(Y, z_a, z_b, null = 1){
  # This function uses glm()-package to estimate the alternative
  # parameters. Right now we need the 4-tuple.
  x <- rep(c(-1,1), length(z_b))
  # This makes it so that z = (z_a,1, z_b,1, z_a,2, z_b,2,....)
  z <- c(t(matrix(c(z_a, z_b), ncol = 2)))
  
  # glm-function to compute the MLE-coefficients. The
  # '-1' in the function removes the standard intercept
  # for the glm-model.
  if(null == 0){ # under the null we do not observe x.
    MLE <- glm(Y ~ z - 1,
               family = 'binomial')$coefficients
  }
  else{
    MLE <- glm(Y ~ x + z - 1,
               family = 'binomial')$coefficients
  }
  
  # note: the coefficients of glm are exactly the 
  # MLE coefficients for the alternative. I also
  # remove the names of the given variables.
  return(unname(MLE / 2))
}


Compute_E_i <- function(z_a_i, 
                        z_b_i, 
                        est_theta_x, 
                        est_theta_1, 
                        Y_i){
  # To compute the E_i-th value, we require:
  # (i) P_1 where theta_x,theta_1 is est. from the (i-1)-data points.
  alternative <- parametrized_alternative(z_a_i, 
                                          z_b_i, 
                                          est_theta_x, 
                                          est_theta_1)
  P_1 <- P_mu(alternative[1], Y_i[1], 
              alternative[2], Y_i[2])
  
  # (ii) The best mixture compared to the est. alternative points.
  band <- P_mixture(z_a_i, 
                    z_b_i, 
                    est_theta_x, # hat{theta}_x^(i-1)
                    est_theta_1, # hat{theta}_1^(i-1)
                    Y_i[1], Y_i[2])
  # The band also holds the P_mix outcome.
  P_mix <- band$P_mix
  
  # (iii): E_i = P_(1,(i-1)) / P_mix.
  return(P_1 / P_mix)
}



