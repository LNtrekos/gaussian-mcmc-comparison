# Needed Libraries:
{
  library(mvtnorm)
  library(coda)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(doSNOW)
}
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
# Parallel Function:
{
  parallel_mcmc <- function(algorithm,
     n_chains = 10, n_cores = 15,
     I = 5e4, burn_in = 1000,
     d = 10, Sigma, gamma = 1/d, thin = 1) {
    
    suppressPackageStartupMessages({
      library(doSNOW)
      library(foreach)
      library(mvtnorm)
    })
    
    # Set up cluster
    cl <- makeCluster(n_cores)
    registerDoSNOW(cl)
    
    # Progress bar setup
    pb <- txtProgressBar(max = n_chains, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    
    # Export variables and load source
    clusterExport(cl, c("d", "Sigma", "gamma", "thin"), envir = environment())
    clusterEvalQ(cl, {
      library(mvtnorm)
      source("C:/Users/lntre/Documents/University/Simulation/Sim_Assigns/Assign_6/mcmc_functions.R")
    })
    
    # Select sampler
    sampler <- switch(algorithm,
      "rwm" = MCMC_sampler_rwm,
      "mala" = MCMC_sampler_mala,
      "gimala" = MCMC_sampler_gimala,
      stop("Invalid algorithm"))
    
    # Run chains
    chains_list <- foreach(
      i = 1:n_chains,
      .packages = "mvtnorm",
      .options.snow = list(progress = progress)
    ) %dopar% {
      
      set.seed(i)
      initial_values <- rnorm(d, mean = 0, sd = 1)
      
      chain <- sampler(dGd_0 = initial_values,
        mcmc_iterations = I,
        Sigma = Sigma, gamma = gamma, thin = thin,
        plot.Trace = FALSE)
      
      keep_samples <- seq(burn_in + 1, I, by = thin)
      chain$dGd_MCMC[keep_samples, ]
    }
    
    # Cleanup
    close(pb)
    stopCluster(cl)  # Explicitly stop cluster
    gc()             # Force garbage collection
    
    # Convert to mcmc.list
    return(mcmc.list(lapply(chains_list, mcmc)))
  }
}
# Example with d = 10:
{
  d = 10
  si = 0.01 * (1:d)
  S = diag(si^2)
  
  # RWM: 
  {
    # Run RWM once : (3.5)
    {
      par(mfrow = c(1,1))
      
      I = 1 * 10^4
      sel.col = round(seq(1, d, length.out = min(5, d)))
      
      rwm_run1 <- MCMC_sampler_rwm(
        dGd_0 = rep(3, d),mcmc_iterations = I,
        gamma = 3.5 * 1/d, Sigma = S,
        thin = 10, plot.Trace = T)
      
      rwm_Burn = 2000
      # Convert to coda object
      rwm_mcmc_chain <- mcmc(rwm_run1$dGd_MCMC[rwm_Burn:I,])
      
      # Trace plot
      matplot(rwm_mcmc_chain[,sel.col], 
              type = 'l', lty = sel.col, col = sel.col,
              xlab = 'Iteration', ylab = 'Parameter Value',
              main = paste('Trace Plot'))
      
      legend("topleft",
             legend = paste("Dim", sel.col),
             col = sel.col, lty = sel.col,
             lwd = 2, cex = 0.8)
      # Density_plots
      par(mfrow = c(2,5))
      densplot(rwm_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Summary statistics
      summary(rwm_mcmc_chain)
      si
      
      # Autocorrelation
      par(mfrow = c(2,5))
      autocorr.plot(rwm_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Effective sample size
      effectiveSize(rwm_mcmc_chain)
      par(mfrow = c(1,1))
    }
    # Multiple RWM Chains:
    {
      # Number of chains
      n_chains <- 10  
      
      # Initialize list to store chains and Iterations:
      chains_list <- vector("list", length = n_chains)
      I = 5 * 10^4
      rwm_Burn = 3000
      
      # Run multiple chains
      for (i in 1:n_chains) {
        
        initial_values <- rnorm(d, mean = 0, sd = 1)  # Random starts
        
        # Run MCMC sampler
        chain <- chain <- MCMC_sampler_rwm(
          dGd_0 = initial_values, mcmc_iterations = I,
          gamma = 1/d, Sigma = S,
          thin = NULL, plot.Trace = FALSE)
        
        # Convert to mcmc object (after burn-in)
        chains_list[[i]] <- mcmc(chain$dGd_MCMC[rwm_Burn:I, ])
        cat(paste0('Chain: ', i, ", Just Finished", '\n'))
      }
      
      multiple_rwm_chains <- mcmc.list(chains_list)
      
      par(mfrow = c(2, 5))
      traceplot(multiple_rwm_chains[, sel.col])
      par(mfrow = c(1,1))
      
      gelman.diag(multiple_rwm_chains)
      
      # Autocorrelation
      autocorr.plot(multiple_rwm_chains)
      
      # Effective sample size
      effectiveSize(multiple_rwm_chains)
      
    }
  }
  # MALA:
  {
    # Run MALA once : (6.9)
    {
      par(mfrow = c(1,1))
      
      I = 1 * 10^4
      sel.col = round(seq(1, d, length.out = min(5, d)))
      
      mala_run1 <- MCMC_sampler_mala(
        dGd_0 = rep(3, d), mcmc_iterations = I,
        gamma = 6.9 *  1/d, Sigma = S, thin = 10, plot.Trace = T)
      
      mala_Burn = 1000
      # Convert to coda object
      mala_mcmc_chain <- mcmc(mala_run1$dGd_MCMC[mala_Burn:I,])
      
      # Trace plot
      matplot(mala_mcmc_chain[, sel.col], 
              type = 'l', lty = sel.col, col = sel.col,
              xlab = 'Iteration', ylab = 'Parameter Value',
              main = paste('Trace Plot'))
      
      legend("topleft",
             legend = paste("Dim", sel.col),
             col = sel.col, lty = sel.col,
             lwd = 2, cex = 0.8)
      
      # Density_plots
      par(mfrow = c(2,5))
      densplot(mala_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Summary statistics
      summary(mala_mcmc_chain)
      
      # Autocorrelation
      par(mfrow = c(2,5))
      autocorr.plot(mala_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Effective sample size
      effectiveSize(mala_mcmc_chain)
      par(mfrow = c(1,1))
    }
    # Multiple MALA Chains:
    {
      # Number of chains
      n_chains <- 10  
      
      # Initialize list to store chains and Iterations:
      chains_list <- vector("list", length = n_chains)
      I = 5 * 10^4
      mala_Burn = 1000
      
      # Run multiple chains
      for (i in 1:n_chains) {
        
        initial_values <- rnorm(d, mean = 0, sd = 1)  # Random starts
        
        # Run MCMC sampler
        chain <- chain <- MCMC_sampler_mala(
          dGd_0 = initial_values, mcmc_iterations = I,
          gamma = 1/d, Sigma = S,
          thin = NULL, plot.Trace = FALSE)
        
        # Convert to mcmc object (after burn-in)
        chains_list[[i]] <- mcmc(chain$dGd_MCMC[mala_Burn:I, ])
        cat(paste0('Chain: ', i, ", Just Finished", '\n'))
      }
      
      multiple_mala_chains <- mcmc.list(chains_list)
      
      par(mfrow = c(2, 5))
      traceplot(multiple_mala_chains[, sel.col])
      par(mfrow = c(1,1))
      
      gelman.diag(multiple_mala_chains)
      
      # Autocorrelation
      autocorr.plot(multiple_mala_chains)
      
      # Effective sample size
      effectiveSize(multiple_mala_chains)
      
    }
  }
  # GI - MALA:
  {
    # Run GI - MALA once :
    {
      par(mfrow = c(1,1))
      
      I = 1 * 10^4
      sel.col = round(seq(1, d, length.out = min(5, d)))
      
      gi_mala_run1 <- MCMC_sampler_gimala(
        dGd_0 = rep(3, d), mcmc_iterations = I,
        gamma = 1, Sigma = S, thin = 1, plot.Trace = T)
      
      gimala_Burn = 1000
      # Convert to coda object
      gimala_mcmc_chain <- mcmc(gi_mala_run1$dGd_MCMC[gimala_Burn:I,])
      
      # Trace plot
      matplot(gimala_mcmc_chain[, sel.col], 
              type = 'l', lty = sel.col, col = sel.col,
              xlab = 'Iteration', ylab = 'Parameter Value',
              main = paste('Trace Plot'))
      
      legend("topleft",
             legend = paste("Dim", sel.col),
             col = sel.col, lty = sel.col,
             lwd = 2, cex = 0.8)
      
      # Density_plots
      par(mfrow = c(2,5))
      densplot(gimala_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Summary statistics
      summary(gimala_mcmc_chain)
      
      # Autocorrelation
      par(mfrow = c(2,5))
      autocorr.plot(gimala_mcmc_chain)
      par(mfrow = c(1,1))
      
      # Effective sample size
      effectiveSize(gimala_mcmc_chain)
      par(mfrow = c(1,1))
    }
    # Multiple GI - MALA Chains:
    {
      # Number of chains
      n_chains <- 10  
      
      # Initialize list to store chains and Iterations:
      chains_list <- vector("list", length = n_chains)
      I = 5 * 10^4
      gimala_Burn = 1000
      
      # Run multiple chains
      for (i in 1:n_chains) {
        
        initial_values <- rnorm(d, mean = 0, sd = 1)  # Random starts
        
        # Run MCMC sampler
        chain <- chain <- MCMC_sampler_gimala(
          dGd_0 = initial_values, mcmc_iterations = I,
          gamma = 1/d, Sigma = S,
          thin = NULL, plot.Trace = FALSE)
        
        # Convert to mcmc object (after burn-in)
        chains_list[[i]] <- mcmc(chain$dGd_MCMC[gimala_Burn:I, ])
        cat(paste0('Chain: ', i, ", Just Finished", '\n'))
      }
      
      multiple_gimala_chains <- mcmc.list(chains_list)
      
      par(mfrow = c(2, 5))
      traceplot(multiple_gimala_chains[, sel.col])
      par(mfrow = c(1,1))
      
      gelman.diag(multiple_gimala_chains)
      
      # Autocorrelation
      autocorr.plot(multiple_gimala_chains)
      
      # Effective sample size
      effectiveSize(multiple_gimala_chains)
      
    }
  }
  # Parallel
  {
    d = 10
    si = 0.01 * (1:d)
    S = diag(si^2)
    
    # Run RWM in parallel
    rwm_Burn = 2000
    rwm_chains <- parallel_mcmc(algorithm = "rwm",
      n_chains = 15, n_cores = 15, I = 5e4, gamma = 3.5 * 1/d, 
      burn_in = rwm_Burn, d = d, Sigma = S)
    
    # Check diagnostics
    summary(rwm_chains)
    gelman.diag(rwm_chains)
    effectiveSize(rwm_chains)
    
    # Run MALA in parallel
    mala_Burn = 1000
    mala_chains <- parallel_mcmc(algorithm = "mala",
      n_chains = 15, n_cores = 15, I = 5e4, gamma = 6.9 * 1/d, 
      burn_in = rwm_Burn, d = d, Sigma = S)
    # Check diagnostics
    summary(mala_chains)
    gelman.diag(mala_chains)
    effectiveSize(mala_chains)
    
    # Run GI-MALA in parallel
    gimala_Burn = 1000
    gimala_chains <- parallel_mcmc(algorithm = "gimala",
       n_chains = 15, n_cores = 15,
       d = d, Sigma = S, burn_in = gimala_Burn)
    
    # Check diagnostics
    summary(gimala_chains)
    gelman.diag(gimala_chains)
    effectiveSize(gimala_chains)
  }
}
# d = 100
# Example with d = 10:
{
  d = 100
  si = 0.01 * (1:d)
  S = diag(si^2)
  
  # RWM (3.5): 
  {
    par(mfrow = c(1,1))
    
    I = 1 * 10^4
    sel.col = round(seq(1, d, length.out = 6))
    
    rwm_run1 <- MCMC_sampler_rwm(
      dGd_0 = rep(3, d),mcmc_iterations = I,
      gamma = 3.5 * 1/d, Sigma = S,
      thin = 10, plot.Trace = T)
    
    rwm_Burn = 5000
    # Convert to coda object
    rwm_mcmc_chain <- mcmc(rwm_run1$dGd_MCMC[rwm_Burn:I,])
    
    # Trace plot
    matplot(rwm_mcmc_chain[,sel.col], 
            type = 'l', lty = sel.col, col = sel.col,
            xlab = 'Iteration', ylab = 'Parameter Value',
            main = paste('Trace Plot'))
    
    legend("topleft",
           legend = paste("Dim", sel.col),
           col = sel.col, lty = sel.col,
           lwd = 2, cex = 0.8)
    # Density_plots
    par(mfrow = c(2,3))
    densplot(rwm_mcmc_chain[, sel.col])
    par(mfrow = c(1,1))
    
    # Summary statistics
    summary(rwm_mcmc_chain)
    
    # Autocorrelation
    par(mfrow = c(2,3))
    autocorr.plot(rwm_mcmc_chain[, sel.col])
    par(mfrow = c(1,1))
    
    # Effective sample size
    effectiveSize(rwm_mcmc_chain)
    par(mfrow = c(1,1))
  } 
  # MALA (30):
  {
    par(mfrow = c(1,1))
    
    I = 1 * 10^4
    sel.col = round(seq(1, d, length.out = 6))
    
    mala_run1 <- MCMC_sampler_mala(
      dGd_0 = rep(3, d), mcmc_iterations = I,
      gamma = 30 *  1/d, Sigma = S, thin = 10, plot.Trace = T)
    
    mala_Burn = 1000
    # Convert to coda object
    mala_mcmc_chain <- mcmc(mala_run1$dGd_MCMC[mala_Burn:I,])
    
    # Trace plot
    matplot(mala_mcmc_chain[, sel.col], 
            type = 'l', lty = sel.col, col = sel.col,
            xlab = 'Iteration', ylab = 'Parameter Value',
            main = paste('Trace Plot'))
    
    legend("topleft",
           legend = paste("Dim", sel.col),
           col = sel.col, lty = sel.col,
           lwd = 2, cex = 0.8)
    
    # Density_plots
    par(mfrow = c(2,3))
    densplot(mala_mcmc_chain[, sel.col])
    par(mfrow = c(1,1))
    
    # Summary statistics
    summary(mala_mcmc_chain)
    
    # Autocorrelation
    par(mfrow = c(2,3))
    autocorr.plot(mala_mcmc_chain[, sel.col])
    par(mfrow = c(1,1))
    
    # Effective sample size
    effectiveSize(mala_mcmc_chain)
    par(mfrow = c(1,1))
  } 
  # GI - MALA (1):
  {
    par(mfrow = c(1,1))
    
    I = 1 * 10^4
    sel.col = round(seq(1, d, length.out = 6))
    
    gi_mala_run1 <- MCMC_sampler_gimala(
      dGd_0 = rep(3, d), mcmc_iterations = I,
      gamma = 1, Sigma = S, thin = 10, plot.Trace = T)
    
    gimala_Burn = 1000
    # Convert to coda object
    gimala_mcmc_chain <- mcmc(gi_mala_run1$dGd_MCMC[gimala_Burn:I,])
    
    # Trace plot
    matplot(gimala_mcmc_chain[, sel.col], 
            type = 'l', lty = sel.col, col = sel.col,
            xlab = 'Iteration', ylab = 'Parameter Value',
            main = paste('Trace Plot'))
    
    legend("topleft",
           legend = paste("Dim", sel.col),
           col = sel.col, lty = sel.col,
           lwd = 2, cex = 0.8)
    
    # Density_plots
    par(mfrow = c(2,3))
    densplot(gimala_mcmc_chain[,sel.col])
    par(mfrow = c(1,1))
    
    # Summary statistics
    summary(gimala_mcmc_chain)
    
    # Autocorrelation
    par(mfrow = c(2, 3))
    autocorr.plot(gimala_mcmc_chain[,sel.col])
    par(mfrow = c(1,1))
    
    # Effective sample size
    effectiveSize(gimala_mcmc_chain)
    par(mfrow = c(1,1))
  }
}
# Parallel
{
  d = 100
  si = 0.01 * (1:d)
  S = diag(si^2)
  
  # Run RWM in parallel
  rwm_Burn = 5000
  rwm_chains <- parallel_mcmc(algorithm = "rwm",
                              n_chains = 15, n_cores = 15, I = 10e4,
                              gamma = 3.5 * 1/d, burn_in = rwm_Burn,
                              d = d, Sigma = S)
  
  # Check diagnostics
  summary(rwm_chains)
  gelman.diag(rwm_chains)
  effectiveSize(rwm_chains)
  
  # Run MALA in parallel
  mala_Burn = 1000
  mala_chains <- parallel_mcmc(algorithm = "mala",
                               n_chains = 15, n_cores = 15, I = 10e4,
                               gamma = 30 * 1/d, burn_in = mala_Burn,
                               d = d, Sigma = S)
  # Check diagnostics
  summary(mala_chains)
  gelman.diag(mala_chains)
  effectiveSize(mala_chains)
  
  # Run GI-MALA in parallel
  gimala_Burn = 1000
  gimala_chains <- parallel_mcmc(algorithm = "gimala",
                                 n_chains = 10, n_cores = 15, I = 10e4,
                                 gamma = 1, burn_in = gimala_Burn,
                                 d = d, Sigma = S)
  
  # Check diagnostics
  summary(gimala_chains)
  gelman.diag(gimala_chains)
  effectiveSize(gimala_chains)
  
  
  library(ggplot2)
  ess_data <- data.frame(
    Algorithm = rep(c("RWM", "MALA", "GI-MALA"), each = 100),
    ESS = c(effectiveSize(rwm_chains), effectiveSize(mala_chains), effectiveSize(gimala_chains))
  )
  ggplot(ess_data, aes(x = Algorithm, y = ESS, fill = Algorithm)) + 
    geom_boxplot() + 
    scale_y_log10() + 
    labs(title = "Effective Sample Size (ESS) in 100 Dimensions", y = "ESS (log scale)") + 
    theme_minimal()
  
  
}
var(sapply(gimala_chains, mean))


N = length(gimala_chains)

