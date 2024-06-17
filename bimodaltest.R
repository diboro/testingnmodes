# ATENTION: this script may last multiple hours to run. Run carefully.

# Script with the simulations of the likelihood ratio test for
#   H0: f has 2 or less modes.
#   H1: f has more than 2 modes.

rm(list = ls())

library(foreach)
library(doParallel)

# We load the models and the functions to compute
# the critical bandwidth and the cross-validation function
source('modcirc.R')
source('models.R')


# Number of simulations
nsim <- 1000
# Number of bootstrap resamples
B <- 500
# Sample size
n <- 100
# n <- 500
# n <- 1000
# The model we generate the data from
model <- "modc1"
# model <- "modc2"
# model <- "modc3"
# model <- "modc4"
# model <- "modc5"
# model <- "modc6"
# model <- "modc7"
# model <- "modc8"
# model <- "modc9"
# model <- "modc10"
# model <- "modc11"
# model <- "modc12"
# model <- "modc13"
# model <- "modc14"
# model <- "modc15"

# Number of cores:
nc <- detectCores() - 2


###### DATA SIMULATION ######

# Seed
set.seed(3141592)

data <- do.call(paste0("r", model), args = list(n = nsim*n))
dim(data) <- c(nsim, n)



# Sample statistics
h2 <- numeric(nsim)
hmax <- numeric(nsim)
h2max <- numeric(nsim)
est <- numeric(nsim)

# Bootstrap statistics
data.boot <- array(dim = c(nsim, B, n))
h2.boot <- matrix(nrow = nsim, ncol = B)
hmax.boot <- matrix(nrow = nsim, ncol = B)
h2max.boot <- matrix(nrow = nsim, ncol = B)
est.boot <- matrix(nrow = nsim, ncol = B)

# p-value
pval <- numeric(nsim)


##### SIMULATION STUDY #####

# ATENTION: this chunk of code may last multiple hours to run.
# Run with caution

# Starting time
tini <- proc.time()
print('Starting simulation study:')
for (i in 1:nsim){
  print(paste0('   Iteration ', i, ':'))
  tini.i <- proc.time()
  # We compute the three bandwidths
  #   - h2
  h2[i] <- bw.crit.circ(data[i,], mod0 = 2)
  #   - hmax
  opt <- optimize(xloglik.circ, c(0, 2*pi),
                  maximum = T, data = data[i, ])
  hmax[i] <- opt$maximum
  #   - h2max
  opt.h2 <- optimize(xloglik.circ, c(h2[i], 2*pi),
                  maximum = T, data = data[i, ])
  h2max[i] <- opt.h2$maximum
  print('      Bandwidths computed.')
  
  # We compute the value of the statistic
  est[i] <- 2 * (opt$objective - opt.h2$objective)
  print('      Statistic computed.')
  
  print('      Starting bootstrap loop:')
  # We parallelize the bootstrap:
  cl <- parallel::makeCluster(nc)
  doParallel::registerDoParallel(cl)
  aux <- foreach(b = 1:B, .combine = 'rbind', 
                 .packages = c('circular', 'denscirc')) %dopar% {
                   # Seed
                   set.seed(i*b)
                   # The bootstrap sample
                   i.boot <- sample(n, replace = T)
                   data.boot <- data[i, i.boot] + 
                     rwrappednormal(n, sd = h2[i])
                   
                   # We compute the three bandwidths
                   #   - h2
                   h2.boot <- bw.crit.circ(data.boot, mod0 = 2)
                   #   - hmax
                   opt.boot <- optimize(xloglik.circ, c(0, 2*pi),
                                        maximum = T, data = data.boot)
                   hmax.boot <- opt.boot$maximum
                   #   - h2max
                   opt.h2.boot <- optimize(xloglik.circ, c(h2.boot, 2*pi),
                                           maximum = T, data = data.boot)
                   h2max.boot <- opt.h2.boot$maximum
                   
                   # We compute the test statistic
                   est.boot <- 2 * (opt.boot$objective - opt.h2.boot$objective)
                   
                   # We return the bootstrap bandwiths, the boostrap statistic
                   # and the bootstrap sample.
                   c(h2.boot, hmax.boot,
                     h2max.boot, est.boot, data.boot)
                 }
  parallel::stopCluster(cl)
  
  # We save all the things that the foarch returns in
  # each numeric vector
  h2.boot[i,] <- aux[, 1]
  hmax.boot[i, ] <- aux[, 2]
  h2max.boot[i, ] <- aux[, 3]
  est.boot[i, ] <- aux[, 4]
  data.boot[i, ,] <- aux[,5:(n+4)]
  print('      Bootstrap loop ended.')
  
  
  # We compute the p-value of the i-th sample
  pval[i] <- mean(est.boot[i,] >= est[i])
  print('      p-value computed.')
  tend.i <- proc.time()
  tempo.i <- tend.i - tini.i
  print(paste0('      This iteration lasted: ',
               format(tempo.i[3], nsmal = 2), ' seg.'))
}
tend <- proc.time()

ctime <- tend-tini
print('Simulation study ended.')

##### SAVING THE RESULTS #####
print('Saving the results...')

# A file with only the p-values:
save(pval, file = paste0('sim_bimodal_', model, '_', n, '.RData'))

# A file with EVERITHING
save(data, h2, hmax, h2max, est, data.boot, h2.boot, hmax.boot,
     h2max.boot, est.boot, pval, ctime,
     file = paste0('sim_bimodal_', model, '_', n, '_verbose.RData'))
print('Results saved.')
