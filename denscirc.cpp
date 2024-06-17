#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;
using namespace std;


/* Function that computes the (circular) KDE
   with wrapped normal kernel in terms of
   its mean resultant length
   Arguments:
   		- data: a numeric vector with the sample
  		- rho: a real number between 0 and 1. 
			It is the mean resultant length
			of the wrapped normal.
		- x: a numeric vector with the points
			where we want to evaluate the KDE.
		- kmax: an integer. The maximum number
			of terms to approximate the series.
			It is equal to 7 by default.
*/
NumericVector denscircrho(NumericVector data,
                          double rho,
                          NumericVector x,
                          int kmax = 7) {
    // The sample size
    int n = data.size();
    // The number of points we want to
    // evaluate the density
    int m = x.size();
    // The density value
    NumericVector dens(m);
    // An auxiliary variable
    double aux;
    // Coefficients
    double coef[kmax];
    
    // In the first iteration we initialize
    // the coefficients
    for(int k = 0; k < kmax; k++){
        aux = 0;
        coef[k] = pow(rho, pow(k+1, 2));
        for(int i = 0; i < n; i++){
            aux += cos((k+1)*(x[0] - data[i]));
        }
        dens[0] += coef[k] * aux;
    }
    dens[0] = (1 + 2 * dens[0]/n)/(2 * M_PI);
    
    // The remaining iterations
    for(int j = 1; j < m; j++){
        for(int k = 0; k < kmax; k++){
            aux = 0;
            for(int i = 0; i < n; i++){
                aux += cos((k+1)*(x[j] - data[i]));
            }
            dens[j] += coef[k] * aux;
        }
        
        dens[j] = (1 + 2 * dens[j]/n)/(2 * M_PI);
    }
    return dens;
}

/* Function that computes the (circular) KDE
   with wrapped normal kernel in terms of
   the bandwidth
   Arguments:
   		- data: a numeric vector with the sample
  		- sigma: a positive real number. 
			It is the standar deviation of the
			normal distribution.
		- x: a numeric vector with the points
			where we want to evaluate the KDE.
		- kmax: an integer. The maximum number
			of terms to approximate the series.
			It is equal to 2 by default.
*/
NumericVector denscircsigma(NumericVector data,
                            double sigma,
                            NumericVector x,
                            int kmax = 2) {
    // The sample size
    int n = data.size();
    // The number of points we want to
    // evaluate the density
    int m = x.size();
    // The density value
    NumericVector dens(m);
    // Auxiliary variables
    double sigma2 = pow(sigma, 2);
    
    // The triple loop
    for (int j = 0; j < m; j++){
        for (int i = 0; i < n; i++){
            for (int k = -kmax; k < kmax; k++){
                dens[j] += exp(- pow(x[j] - data[i] + 2*k*M_PI, 2)/(2*sigma2));
            }
        }
        dens[j] = dens[j]/(sqrt(2*M_PI*sigma2)*n);
    }
    return dens;
}

/* Function that computes the (circular) KDE
   with wrapped normal kernel in terms of
   the bandwidth
   Arguments:
   		- data: a numeric vector with the sample
  		- bw: the bandwidth. It must be
			a positive real number.
		- x: a numeric vector with the points
			where we want to evaluate the KDE.
*/

// [[Rcpp::export]]
NumericVector denscirc(NumericVector data,
                       double bw,
                       NumericVector x) {
    // The number of points we want to
    // evaluate the density
    int m = x.size();
    // The density value
    NumericVector dens(m);
    
    // We compute the density differently
    // depending on the value of bw
    if (bw <= 1) {
        dens = denscircsigma(data, bw, x);
    } else {
        double rho = exp(-pow(bw, 2)/ 2);
        int kmax = ceil(sqrt(log(1e-12)/log(rho)));
        dens = denscircrho(data, rho, x, kmax);
    }
    
    return dens;
}

/* A function that creates a sequence of evenly-spaced points.
   Arguments:
   		- a, b: two real numbers. They are the starting
			and end value of the sequence, respectively.
  		- N: the number of elements of the sequence.
*/
// [[Rcpp::export]]
NumericVector range(double a, double b, std::size_t N)
{
    double h = (b - a) / static_cast<double>(N-1);
    NumericVector xs(N);
    NumericVector::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

/* A function that computes the number of modes of
   the (circular) KDE for a given bandwidth.
   Arguments:
		- data: a numeric vector with the data sample.
		- bw: a positive real number. The bandwidth.
		- n: an integer. The number of points where 
			we want to evaluate the KDE. 512 by default.
		- dens: a numeric vector of positive real numbers.
			If nodens = false, this vector is used as
			the KDE evaluated in a sorted grid of points.
			If nodens = true, it is ignored.
		- nodens: boolean. It indicates if a vector with
			the KDE evaluated in a sorted grid of points
			is provided or not.
*/

// [[Rcpp::export]]
int nmodcirc(NumericVector data,
             double bw,
             int n = 512,
             NumericVector dens = 0,
             bool nodens = true){
    // Output variable: number of modes
    int nmod = 0;
    
    // First, we check if the density vector is provided
    if (nodens){
        // If the density vector is not provided,
        // we have to compute it
        NumericVector x(n);
        x = range(0, 2*M_PI, n);
        dens = denscirc(data, bw, x);
    } else {
        // If it is provided, we update n so it is equal to its size
        n = dens.size();
    }
    
    // We compute the modes
    for (int i = 1; i < (n-1); i++){
        // If dens[i] is greater than dens[i-1] and dens[i+1]
        // there is a mode
        if ((dens[i] - dens[i-1] > 0) && (dens[i] - dens[i+1] > 0)) {
            nmod++;
        }
    }
    
    // Finaly, we check if there is a mode in dens[0]
    if ((dens[n-1] - dens[n-2] > 0) && (dens[0] - dens[1] > 0)) {
        nmod++;
    }
    return nmod;
}
