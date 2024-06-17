library(circular)
# We need to load all the Rcpp functions in a package
# If not, it is not possible to parallelise the code
# To create a package with all the Rcpp functions
# in "denscirc.cpp" and install it in your computer,
# uncomment and run the following three lines:

# library(Rcpp)
# Rcpp.package.skeleton("denscirc",cpp_files = "denscirc.cpp")
# install.packages("denscirc", repos=NULL, type="source")

# Once the package is installed, it can be loaded as usual package
library(denscirc)

# Function that computes the critical bandwidth for k modes of a circular sample
# The critical bandwidth is computed by a bisection method
# The arguments are:
#   - data: the data sample
#   - mod0: the number of modes k
#   - n: the number of points where the density function is evaluated
#   - bw.ini: the initial iteration of the bandwidth
#   - bw.down, bw.up: the lower and upper extremes of the interval where
#       the critical bandwidth is searched. If they are not provided,
#       they are computed multiplying bw.ini by 0.5 (respectively, by 2)
#       until the number of modes is <= mod0 (respectively, > mod0)
#   - tol: the tolerance. The algorithm stops when (bw.up - bw.down < tol)
bw.crit.circ <- function(data, mod0 = 1, n = 2^9,
                         bw.ini = 1, bw.down = NaN, bw.up = NaN, 
                         tol = 1e-5){
    # If bw.down and bw.up are not provided, we need to compute them
    while(any(is.na(c(bw.down, bw.up)))){
        # Depending on the number of modes of bw.ini, we do different things
        n.ini <- nmodcirc(data, bw.ini)
        if (n.ini > mod0){
            # If the number of modes is greater than mod0,
            # then bw.ini is smaller than the critical bandwidth.
            # We make bw.ini bigger by multiplying it by 2
            bw.down <- bw.ini
            bw.ini <- 2*bw.ini
        } else {
            # If the number of modes is smaller than mod0,
            # then bw.ini is greater than the critical bandwidth.
            # We make bw.ini smaller by multiplying it by 0.5
            bw.up <- bw.ini
            bw.ini <- 0.5*bw.ini
        }
    }
    
    
    # Once bw.down and bw.up are computed, we compute the critical bandwidth
    # by a bisection method. The process stops when (bw.up - bw.down < tol)
    while(bw.up - bw.down > tol){
        # The new bandwidth and its number of modes
        bw.new <- (bw.up + bw.down)/2
        n.new <- nmodcirc(data, bw.new)
        if (n.new <= mod0){
            # If the number of modes is <= mod0, then bw.new is
            # larger than the critical bandwidth. We actualize bw.up
            bw.up <- bw.new
        } else {
            # If the number of modes is > mod0, then bw.new is
            # smaller than the critical bandwidth. We actualize bw.down
            bw.down <- bw.new
        }
    }
    # We return bw.up as an approximation of the bandwidth
    return(bw.up)
}


# Auxiliary function that computes the circular KDE
# without the i-th point of the sample and evaluate it in it
# The arguments are:
#   - i: the position of the point that we want to take
#       out of the sample
#   - h: the bandwidth. Must be a real number.
#   - data: the sample
xfi.circ <- function(i, h, data){
    f <- denscirc(data[-i], h, data[i])
}

# Function that computes the cross-validation (circular) KDE
# in all points of the sample
# The arguments are:
#   - h: the bandwidth. Must be a real number.
#   - data: the sample
xf.circ <- function(h, data){
    # Sample size
    n <- length(data)
    # The kde withouth cross-validation
    f <- denscirc(data, h, data)
    # The value of K_h(0)/n
    pen <- rep(dwrappednormal(circular(0))/(n * h), n)
    
    # ATENTION! The difference (f - pen) may cause underflow problems 
    # when h is small. In that case, we use the function xfi.circ()
    # to evaluate the cross-validation circular KDE
    res <- (n/(n - 1))*(f - pen)
    if (all(res > 0)){
        return(res)
    } else {
        return(sapply(1:n, xfi.circ, h = h, data = data))
    }
}

# Function that computes the cross-validation likelihood function
# as a function of the bandwidth h
# The arguments are:
#   - h: the bandwidth. Must be a real number.
#   - data: the sample
xloglik1.circ <- function(h, data){
    sum(log(xf.circ(h, data)))
}

# Vectorized version of the function xloglik1.circ()
# The arguments are:
#   - h: the bandwidth. It can be a vector of real numbers.
#   - data: the sample
xloglik.circ <- function(h, data){
    sapply(h, xloglik1.circ, data = data)
}