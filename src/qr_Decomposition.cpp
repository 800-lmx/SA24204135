#include <Rcpp.h>
using namespace Rcpp;
//' @name qr_decomposition
//' @title QR Decomposition
//' @description Compute the QR decomposition of a matrix A
//' @param A A numeric matrix
//' @return A list containing matrices Q and R
//' @export
// [[Rcpp::export]]
 List qr_decomposition(NumericMatrix A) {
   int m = A.nrow();
   int n = A.ncol();
   
   NumericMatrix Q(m, n);
   NumericMatrix R(n, n);
   
   // Perform Gram-Schmidt process
   for (int j = 0; j < n; j++) {
     // Set R(j, j) = ||A(:, j)||
     R(j, j) = sqrt(sum(A(_, j) * A(_, j)));
     
     // Set Q(:, j) = A(:, j) / R(j, j)
     Q(_, j) = A(_, j) / R(j, j);
     
     // Orthogonalize the remaining columns
     for (int k = j + 1; k < n; k++) {
       R(j, k) = sum(Q(_, j) * A(_, k));
       A(_, k) = A(_, k) - Q(_, j) * R(j, k);
     }
   }
   
   return List::create(Named("Q") = Q, Named("R") = R);
 }










//' @name monte_carlo_pi
//' @title Monte Carlo Estimation of Pi
//' @description Use Monte Carlo simulation to estimate the value of Pi
//' @param n An integer specifying the number of random points to generate
//' @return An estimate of Pi
//' @export
// [[Rcpp::export]]

 double monte_carlo_pi(int n) {
   int inside_circle = 0;
   
   // Generate n random points
   for (int i = 0; i < n; i++) {
     // Generate random x and y between -1 and 1
     double x = R::runif(-1, 1);
     double y = R::runif(-1, 1);
     
     // Check if the point is inside the unit circle
     if (x * x + y * y <= 1) {
       inside_circle++;
     }
   }
   
   // Estimate Pi
   double pi_estimate = 4.0 * inside_circle / n;
   return pi_estimate;
 }
