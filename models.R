###### Considered Models in the Simulation Study #####

# First, we need to load the package "circular"
library(circular)

#### Unimodal ####

### Model 1: von-Mises vM(pi, 1)
dmodc1 <- function(x) dvonmises(x, circular(pi), 1)
rmodc1 <- function(n) rvonmises(n, circular(pi), 1)

### Model 2: a mixture 0.2vM(2π/3, 3) + 0.6vM(π, 1.4) + 0.2vM(4π/3, 3).
dmodc2 <- function(x){
    von1 <- dvonmises(x, circular(2*pi/3), 3)
    von2 <- dvonmises(x, circular(pi), 1.4)
    von3 <- dvonmises(x, circular(4*pi/3), 3)
    return(0.2*von1 + 0.6*von2 + 0.2*von3)
}
rmodc2 <- function(n) {
    von1 <- rvonmises(n, circular(2*pi/3), 3)
    von2 <- rvonmises(n, circular(pi), 1.4)
    von3 <- rvonmises(n, circular(4*pi/3), 3)
    iwhich <- sample(1:3, n, prob = c(0.2, 0.6, 0.2), replace = T)
    return(von1*(iwhich == 1) + von2*(iwhich == 2) + von3*(iwhich == 3))
}

### Model 3: a mixture 0.05·vM(2π/3, 7) + 0.9·vM(π, 1) + 0.05·vM(4π/3, 7).
dmodc3 <- function(x){
    von1 <- dvonmises(x, circular(2*pi/3), 7)
    von2 <- dvonmises(x, circular(pi), 1)
    von3 <- dvonmises(x, circular(4*pi/3), 7)
    return(0.05*von1 + 0.9*von2 + 0.05*von3)
}
rmodc3 <- function(n) {
    von1 <- rvonmises(n, circular(2*pi/3), 7)
    von2 <- rvonmises(n, circular(pi), 1)
    von3 <- rvonmises(n, circular(4*pi/3), 7)
    iwhich <- sample(1:3, n, prob = c(0.05, 0.9, 0.05), replace = T)
    return(von1*(iwhich == 1) + von2*(iwhich == 2) + von3*(iwhich == 3))
}

### Model 4: a sine-skewed von–Mises: ssvM(π, 1,−0.9).
# Density function of kssvm
dkssvm <- function(x, mu, kappa, lambda, k){
    dvm <- dvonmises(x, mu = mu, kappa = kappa)
    skew <- 1 - lambda * sin(k * (x - mu))
    return(as.numeric(dvm * skew))
}
# Density function of the model
dmodc4 <- function(x){
    dkssvm(x, circular(pi), 1, 0.9, 1)
}
# Random function of kssvm
rkssvm <- function(n, mu, kappa, lambda, k){
    x <- rvonmises(n, mu = mu, kappa = kappa)
    u <- runif(n)
    neg <- dkssvm(x, mu, kappa, lambda, k) < u * dvonmises(x, mu, kappa)
    return(neg * (2*pi - x) + (1 - neg) * x)
}
# Random function of the model
rmodc4 <- function(n) rkssvm(n, circular(pi), 1, 0.9, 1)

### Model 5: a wrapped beta wBeta(3,2)
# Density function of wBeta(shape1, shape2) 
# with support (pi/2, 3pi/2)
dwbeta <- function(x, shape1, shape2){
    s <- (as.numeric(x) - 0.5*pi)/pi
    dbeta(s, shape1, shape2)/pi
}
# Density function of the model
dmodc5 <- function(x){
    dwbeta(x, 3, 2)
}
# Random function of the model
rmodc5 <- function(n){
    circular(rbeta(n, 3, 2)*pi + pi/2)
}

# All five together:
plotunimod <- function(){
    par(mfrow=c(2,3))
    densities <- list(dmodc1, dmodc2, dmodc3, dmodc4, dmodc5)
    for(i in 1:5){
        curve(densities[[i]](circular(x)), from = 0, to = 2*pi,
              ylab = paste('Circular Model', i))
    } 
}

#### Bimodal ####

### Model 6: a mixture of two von–Mises: 0.5·vM(2, 3) + 0.5·vM(4, 3).
dmodc6 <- function(x){
    von1 <- dvonmises(x, circular(pi-1.25), 1.5)
    von2 <- dvonmises(x, circular(pi+1.25), 1.5)
    return(0.5*von1 + 0.5*von2)
}
rmodc6 <- function(n){
    von <- rvonmises(n, circular(pi), 1.5)
    pos <- sample(c(-1,1), n, replace = T)
    return(von + 1.25*pos)
}

### Model 7: a mixture of von–Mises: 0.5·vM(π − 1, 1.5) + 0.5·vM(π + 1, 1.5).
dmodc7<- function(x){
    von1 <- dvonmises(x, circular(pi - 1), 1.5)
    von2 <- dvonmises(x, circular(pi + 1), 1.5)
    return(0.5*von1 + 0.5*von2)
}
rmodc7 <- function(n){
    von <- rvonmises(n, circular(pi), 1.5)
    pos <- sample(c(-1,1), n, replace = T)
    return(von + pos)
}

### Model 8: a mixture of two von–Mises:
### 0.5·vM(1.5, 4) + 0.5·vM(3, 2).
dmodc8 <- function(x){
    von1 <- dvonmises(x, circular(1.5), 4)
    von2 <- dvonmises(x, circular(3), 2)
    return(0.5*von1 + 0.5*von2)
}
rmodc8 <- function(n){
    von1 <- rvonmises(n, circular(1.5), 4)
    von2 <- rvonmises(n, circular(3), 2)
    iwhich <- sample(0:1, n, replace = T)
    return(von1*iwhich + von2*(1-iwhich))
}

### Model 9: a mixture of two von–Mises: 0.95·vM(π/2, 2) + 0.05·vM(3π/2, 3).
dmodc9 <- function(x){
    von1 <- dvonmises(x, circular(pi/2), 2)
    von2 <- dvonmises(x, circular(3*pi/2), 3)
    return(0.95*von1 + 0.05*von2)
}
rmodc9 <- function(n){
    von1 <- rvonmises(n, circular(pi/2), 2)
    von2 <- rvonmises(n, circular(3*pi/2), 3)
    iwhich <- sample(0:1, n, replace = T, prob = c(0.05, 0.95))
    return(von1*iwhich + von2*(1-iwhich))
}

### Model 10: a mixture of two von–Mises: 0.9·vM(π/2, 2) + 0.10·vM(3π/2, 3).
dmodc10 <- function(x){
    von1 <- dvonmises(x, circular(pi/2), 2)
    von2 <- dvonmises(x, circular(3*pi/2), 3)
    return(0.9*von1 + 0.1*von2)
}
rmodc10 <- function(n){
    von1 <- rvonmises(n, circular(pi/2), 2)
    von2 <- rvonmises(n, circular(3*pi/2), 3)
    iwhich <- sample(0:1, n, replace = T, prob = c(0.1, 0.9))
    return(von1*iwhich + von2*(1-iwhich))
}

# All four together:
plotbimod <- function(){
    par(mfrow=c(2,3))
    densities <- list(dmodc6, dmodc7, dmodc8, dmodc9, dmodc10)
    for(i in 1:5){
        curve(densities[[i]](circular(x)), from = 0, to = 2*pi,
              ylab = paste('Circular Model', i+5))
    }
}

#### Trimodal ####

### Model 11: a mixture of three von–Mises: 
# 1/3·vM(pi-2, 7) + 1/3·vM(pi, 7) + 1/3·VM(pi+2, 7)
dmodc11 <- function(x){
    von1 <- dvonmises(x, circular(pi-2), 7)
    von2 <- dvonmises(x, circular(pi), 7)
    von3 <- dvonmises(x, circular(pi+2), 7)
    return((von1 + von2 + von3)/3)
}

rmodc11 <- function(n){
    von1 <- rvonmises(n, circular(pi-2), 7)
    von2 <- rvonmises(n, circular(pi), 7)
    von3 <- rvonmises(n, circular(pi+2), 7)
    pos <- sample(1:3, n, replace = T, prob = rep(1/3, 3))
    res <- ifelse(pos == 1, von1, ifelse(pos == 2, von2, von3))
    return(res)
}

### Model 12: a mixture of four von–Mises: 
# 1/3·vM(pi-1, 7) + 1/3·vM(pi, 7) + 1/3·VM(pi+1, 7)
dmodc12 <- function(x){
    von1 <- dvonmises(x, circular(pi-1), 7)
    von2 <- dvonmises(x, circular(pi), 7)
    von3 <- dvonmises(x, circular(pi+1), 7)
    return((von1 + von2 + von3)/3)
}

rmodc12 <- function(n){
    von1 <- rvonmises(n, circular(pi-1), 7)
    von2 <- rvonmises(n, circular(pi), 7)
    von3 <- rvonmises(n, circular(pi+1), 7)
    pos <- sample(1:3, n, replace = T, prob = rep(1/3, 3))
    res <- ifelse(pos == 1, von1, ifelse(pos == 2, von2, von3))
    return(res)
}

### Model 13: a mixture of three von–Mises: 
# 0.2·vM(pi/2, 6) + 0.2·vM(7pi/8, 6) + 0.6·vM(7pi/4, 8)
dmodc13 <- function(x){
    von1 <- dvonmises(x, circular(pi/2), 6)
    von2 <- dvonmises(x, circular(pi), 6)
    von3 <- dvonmises(x, circular(7*pi/4), 8)
    return(0.2*von1 + 0.2*von2 + 0.4*von3)
}

rmodc13 <- function(n){
    von1 <- rvonmises(n, circular(pi/2), 6)
    von2 <- rvonmises(n, circular(pi), 6)
    von3 <- rvonmises(n, circular(7*pi/4), 8)
    pos <- sample(1:3, n, replace = T, prob = c(0.2, 0.2, 0.4))
    res <- ifelse(pos == 1, von1, ifelse(pos == 2, von2, von3))
    return(res)
}

### Model 14: a mixture of three von–Mises: 
# 0.1·vM(pi/2, 6) + 0.25·vM(pi, 6) + 0.65·vM(7pi/4, 8).
dmodc14 <- function(x){
    von1 <- dvonmises(x, circular(pi/2), 6)
    von2 <- dvonmises(x, circular(pi), 6)
    von3 <- dvonmises(x, circular(7*pi/4), 8)
    return(0.1*von1 + 0.25*von2 + 0.45*von3)
}

rmodc14 <- function(n){
    von1 <- rvonmises(n, circular(pi/2), 6)
    von2 <- rvonmises(n, circular(pi), 6)
    von3 <- rvonmises(n, circular(7*pi/4), 8)
    pos <- sample(1:3, n, replace = T, prob = c(0.1, 0.25, 0.45))
    res <- ifelse(pos == 1, von1, ifelse(pos == 2, von2, von3))
    return(res)
}

### Model 15: a mixture of three von–Mises: 
# 0.2·vM(pi/2, 6) + 0.2·vM(6pi/7, 6) + 0.6·vM(7pi/4, 8)
dmodc15 <- function(x){
    von1 <- dvonmises(x, circular(pi/2), 6)
    von2 <- dvonmises(x, circular(6*pi/7), 6)
    von3 <- dvonmises(x, circular(7*pi/4), 8)
    return(0.2*von1 + 0.2*von2 + 0.4*von3)
}

rmodc15 <- function(n){
    von1 <- rvonmises(n, circular(pi/2), 6)
    von2 <- rvonmises(n, circular(6*pi/7), 6)
    von3 <- rvonmises(n, circular(7*pi/4), 8)
    pos <- sample(1:3, n, replace = T, prob = c(0.2, 0.2, 0.4))
    res <- ifelse(pos == 1, von1, ifelse(pos == 2, von2, von3))
    return(res)
}

# All four together:
plottrimod <- function(){
    par(mfrow=c(2,3))
    densities <- list(dmodc11, dmodc12, dmodc13, dmodc14, dmodc15)
    for(i in 1:5){
        curve(densities[[i]](circular(x)), from = 0, to = 2*pi,
              ylab = paste('Circular Model', i+10))
    }
}
