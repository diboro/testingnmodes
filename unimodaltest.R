# ATENTION: this script may last multiple hours to run. Run carefully.

# Script with the simulations of the likelihood ratio test for
#   H0: f has 1 mode.
#   H1: f has more than 1 mode.
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
h1 <- numeric(nsim)
hmax <- numeric(nsim)
h1max <- numeric(nsim)
est <- numeric(nsim)

# Bootstrap statistics
data.boot <- array(dim = c(nsim, B, n))
h1.boot <- matrix(nrow = nsim, ncol = B)
hmax.boot <- matrix(nrow = nsim, ncol = B)
h1max.boot <- matrix(nrow = nsim, ncol = B)
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
  #   - h1
  h1[i] <- bw.crit.circ(data[i,], mod0 = 1)
  #   - hmax
  opt <- optimize(xloglik.circ, c(0, 2*pi),
                  maximum = T, data = data[i, ])
  hmax[i] <- opt$maximum
  #   - h1max
  opt.h1 <- optimize(xloglik.circ, c(h1[i], 2*pi),
                  maximum = T, data = data[i, ])
  h1max[i] <- opt.h1$maximum
  print('      Bandwidths computed.')
  
  # We compute the value of the statistic
  est[i] <- 2 * (opt$objective - opt.h1$objective)
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
                     rwrappednormal(n, sd = h1[i])
                   
                   # We compute the three bandwidths
                   #   - h1
                   h1.boot <- bw.crit.circ(data.boot)
                   #   - hmax
                   opt.boot <- optimize(xloglik.circ, c(0, 2*pi),
                                        maximum = T, data = data.boot)
                   hmax.boot <- opt.boot$maximum
                   #   - h1max
                   opt.h1.boot <- optimize(xloglik.circ, c(h1.boot, 2*pi),
                                           maximum = T, data = data.boot)
                   h1max.boot <- opt.h1.boot$maximum
                   
                   # We compute the test statistic
                   est.boot <- 2 * (opt.boot$objective - opt.h1.boot$objective)
                   
                   # We return the bootstrap bandwiths, the boostrap statistic
                   # and the bootstrap sample.
                   c(h1.boot, hmax.boot,
                     h1max.boot, est.boot, data.boot)
                 }
  parallel::stopCluster(cl)
  
  # We save all the things that the foarch returns in
  # each numeric vector
  h1.boot[i,] <- aux[, 1]
  hmax.boot[i, ] <- aux[, 2]
  h1max.boot[i, ] <- aux[, 3]
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
save(pval, file = paste0('sim_', model, '_', n, '.RData'))

# A file with EVERITHING
save(data, h1, hmax, h1max, est, data.boot, h1.boot, hmax.boot,
     h1max.boot, est.boot, pval, ctime,
     file = paste0('sim_', model, '_', n, '_verbose.RData'))
print('Results saved.')
