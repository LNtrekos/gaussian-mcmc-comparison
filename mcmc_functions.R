# Algorithms:
{
  ######
  # RWM: 
  MCMC_sampler_rwm <- function(
    dGd_0 = rnorm(d, mean = 0, sd = 1),
    mcmc_iterations,
    Sigma = NULL, gamma = 1/d,
    thin = 10, plot.Trace = TRUE){
    
    # Initialazation Step:
    d <- length(dGd_0)
    dGdSim <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    ergodic_mean <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    acceptance_rate <- 0
    
    # Precompute precision matrix
    Sigma_inv <- solve(Sigma)
    
    # Initial state
    dGdSim[1, ] <- dGd_0
    ergodic_mean[1, ] <- dGd_0
    
    for (i in 2:mcmc_iterations) {
      # Propose new state
      dGdProposed <- t(rmvnorm(1, mean = dGdSim[i - 1, ],
                               sigma = 2 * gamma * Sigma))
      
      # Log-acceptance ratio
      log_acceptance_ratio <- -0.5 * (
        as.numeric(t(dGdProposed) %*% Sigma_inv %*% dGdProposed) -
          as.numeric(t(dGdSim[i-1, ]) %*% Sigma_inv %*% dGdSim[i-1, ]))
      log_acceptance_ratio <- min(0, log_acceptance_ratio)
      
      # Metropolis step
      if (log(runif(1)) < log_acceptance_ratio) {
        dGdSim[i, ] <- dGdProposed
        acceptance_rate <- acceptance_rate + 1
      } else {
        dGdSim[i, ] <- dGdSim[i-1, ]
      }
      
      # Ergodic mean update
      ergodic_mean[i, ] <- ((i - 1) * ergodic_mean[i - 1, ] + dGdSim[i, ]) / i
      
      # Plotting (subsampled)
      if (plot.Trace && i %% 100 == 0) {
        cat(paste0('iteration: ', i, ', a.r = ', 
                   round(100 * acceptance_rate/i, 2), "%"), '\n')
        
        index_subset <- seq(1, i, by = thin)
        
        # Select 5 dimensions to plot (for clarity)
        dims_to_plot <- round(seq(1, d, length.out = min(5, d)))
        
        # Trace plot
        matplot(index_subset, 
                dGdSim[index_subset, dims_to_plot], 
                type = 'l', lty = 1:d, col = 1:d,
                xlab = 'Iteration', ylab = 'Parameter Value',
                main = paste('Trace Plot (Thinned: every', thin, 'iterations)'))
        
        legend("topleft",
               legend = paste("Dim", dims_to_plot),
               col = 1:d,lty = 1:d, lwd = 2, cex = 0.8)
      }
    }
    
    return(list(
      dGd_MCMC = dGdSim,
      accept_rate = acceptance_rate / mcmc_iterations,
      theta_ergodic_mean = ergodic_mean
    ))
  }
  #######
  # MALA:
  MCMC_sampler_mala <- function(
    dGd_0 = rnorm(d, mean = 0, sd = 1),
    mcmc_iterations,
    Sigma = NULL, gamma = 1/d,
    thin = 10, plot.Trace = TRUE){
    
    # Initialization Step:
    d <- length(dGd_0)
    dGdSim <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    ergodic_mean <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    acceptance_rate <- 0
    
    # Precompute precision matrix
    Sigma_inv <- solve(Sigma)
    
    # Initial state
    dGdSim[1, ] <- dGd_0
    ergodic_mean[1, ] <- dGd_0
    
    for (i in 2:mcmc_iterations) {
      # Propose new state
      dGdProposed <- t(rmvnorm(1, mean = (1 - gamma) * dGdSim[i - 1, ],
                               sigma = 2 * gamma * Sigma))
      
      # Log-acceptance ratio
      log_acceptance_ratio <- - (gamma/4) * (
        as.numeric(t(dGdProposed) %*% Sigma_inv %*% dGdProposed) -
          as.numeric(t(dGdSim[i-1, ]) %*% Sigma_inv %*% dGdSim[i-1, ]))
      
      log_acceptance_ratio <- min(0, log_acceptance_ratio)
      
      # Metropolis step
      if (log(runif(1)) < log_acceptance_ratio) {
        dGdSim[i, ] <- dGdProposed
        acceptance_rate <- acceptance_rate + 1
      } else {
        dGdSim[i, ] <- dGdSim[i-1, ]
      }
      
      # Ergodic mean update
      ergodic_mean[i, ] <- ((i - 1) * ergodic_mean[i - 1, ] + dGdSim[i, ]) / i
      
      # Plotting (subsampled)
      if (plot.Trace && i %% 100 == 0) {
        cat(paste0('iteration: ', i, ', a.r = ', 
                   round(100 * acceptance_rate/i, 2), "%"), '\n')
        
        index_subset <- seq(1, i, by = thin)
        
        # Select 5 dimensions to plot (for clarity)
        dims_to_plot <- round(seq(1, d, length.out = min(5, d)))
        
        # Trace plot
        matplot(index_subset, 
                dGdSim[index_subset, dims_to_plot], 
                type = 'l', lty = 1:d, col = 1:d,
                xlab = 'Iteration', ylab = 'Parameter Value',
                main = paste('Trace Plot (Thinned: every', thin, 'iterations)'))
        
        legend("topleft",
               legend = paste("Dim", dims_to_plot),
               col = 1:d,lty = 1:d, lwd = 2, cex = 0.8)
      }
    }
    
    return(list(
      dGd_MCMC = dGdSim,
      accept_rate = acceptance_rate / mcmc_iterations,
      theta_ergodic_mean = ergodic_mean
    ))
  }
  ############
  # GI - MALA:
  MCMC_sampler_gimala <- function(
    dGd_0 = rnorm(d, mean = 0, sd = 1),
    mcmc_iterations,
    Sigma = NULL, gamma = 1/d,
    thin = 10, plot.Trace = TRUE){
    
    # Initialization Step:
    d <- length(dGd_0)
    dGdSim <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    ergodic_mean <- matrix(NA, nrow = mcmc_iterations, ncol = d)
    acceptance_rate <- 0
    
    # Precompute precision matrix
    Sigma_inv <- solve(Sigma)
    
    # Initial state
    dGdSim[1, ] <- dGd_0
    ergodic_mean[1, ] <- dGd_0
    
    for (i in 2:mcmc_iterations) {
      
      # Propose new state and set:
      dGdSim[i, ] <- t(rmvnorm(1, mean = (1 - gamma) * dGdSim[i - 1, ],
                               sigma = (2*gamma - gamma^2) * Sigma))
      
      # Plotting (subsampled)
      if (plot.Trace && i %% 100 == 0) {
        cat(paste0('iteration: ', i, '\n'))
        
        index_subset <- seq(1, i, by = thin)
        
        # Select 5 dimensions to plot (for clarity)
        dims_to_plot <- round(seq(1, d, length.out = min(5, d)))
        
        # Trace plot
        matplot(index_subset, 
                dGdSim[index_subset, dims_to_plot], 
                type = 'l', lty = 1:d, col = 1:d,
                xlab = 'Iteration', ylab = 'Parameter Value',
                main = paste('Trace Plot (Thinned: every', thin, 'iterations)'))
        
        legend("topleft",
               legend = paste("Dim", dims_to_plot),
               col = 1:d,lty = 1:d, lwd = 2, cex = 0.8)
      }
    }
    
    return(list(
      dGd_MCMC = dGdSim,
      accept_rate = acceptance_rate / mcmc_iterations,
      theta_ergodic_mean = ergodic_mean
    ))
  }
  
}