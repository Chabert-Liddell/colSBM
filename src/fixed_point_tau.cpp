#include <RcppArmadillo.h>
using namespace Rcpp;

// // Function to calculate xlogy
// inline NumericMatrix xlogy(const NumericMatrix &x, const NumericMatrix &y, double eps = 1e-9)
// {
//     return log(x + eps) * y;
// }

// Function to calculate logit
inline arma::mat logit(const arma::mat &x)
{
    return log(x/ (1.0 - x));
}

// // Function to calculate log
// inline NumericMatrix log(const NumericMatrix &x, double eps = 1e-9)
// {
//     return log(x + eps);
// }

// // Function to calculate softmax of a matrix
// NumericMatrix softmax(const NumericMatrix &matrix)
// {
//     int rows = matrix.nrow();
//     int cols = matrix.ncol();
//     NumericMatrix result(rows, cols);

//     for (int i = 0; i < rows; i++)
//     {
//         double maxVal = Rcpp::max(matrix(i, _));
//         double sumExp = 0.0;
//         for (int j = 0; j < cols; j++)
//         {
//             result(i, j) = exp(matrix(i, j) - maxVal);
//             sumExp += result(i, j);
//         }
//         for (int j = 0; j < cols; j++)
//         {
//             result(i, j) /= sumExp;
//         }
//     }

//     return result;
// }

// The fixed point algorithm to update the tau
//          tau_new <-
// if (d == 1)
// {
// #n[[1]] * Q1
//     tau_new < -t(matrix(.xlogy(self$Cpi[[1]][, m],
//                                self$pi [[m]] [[d]], eps = NULL),
//                         self$Q[d], self$n[[1]][m])) +
//                   ((self$nonNAs [[m]]) * self$A [[m]]) % * %
//                       t(self$Cpi[[2]][, m] * t(self$tau [[m]][[2]])) % * %
//                       t(.logit(self$Calpha * self$delta[m] * self$alpha,
//                                eps = 1e-9)) +
//                   (self$nonNAs [[m]]) % * %
//                       t(self$Cpi[[2]][, m] * t(self$tau [[m]][[2]])) % * %
//                       t(.log(1 - self$Calpha * self$alpha * self$delta[m],
//                              eps = 1e-9))
// #In order to fix NaN appearing in the formula(log(Pi) when Pi
// #= 0), the.xlogy function is used with eps = 1e-9
// #POSSIBLE POINT OF FAILURE
// }
// if (d == 2)
// {
// #n[[2]] * Q2
//     tau_new <- t(matrix(.xlogy(self$Cpi[[2]][, m],
//                                self$pi [[m]] [[d]], eps = NULL),
//                         self$Q[d], self$n[[2]][m])) +
//                   t((self$nonNAs [[m]]) * self$A [[m]]) % * %
//                       t(self$Cpi[[1]][, m] * t(self$tau [[m]][[1]])) % * %
//                                                                            .logit(self$Calpha * self$delta[m] * self$alpha, eps = 1e-9) +
//                   t(self$nonNAs [[m]]) % * %
//                       t(self$Cpi[[1]][, m] * t(self$tau [[m]][[1]])) % * %
//                                                                            .log(1 - self$Calpha * self$alpha * self$delta[m], eps = 1e-9)
// #In order to fix NaN appearing in the formula(log(Pi) when Pi
// #= 0), the.xlogy function is used with eps = 1e-9
// #POSSIBLE POINT OF FAILURE
// }
// invisible(tau_new)
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat fixed_point_tau(
    int d, arma::mat tau_m_old_other_dim,
    arma::vec Cpi_1_m, arma::vec Cpi_2_m, 
    arma::mat pi_m_d, int Q_d, int n_d_m,
    arma::mat nonNAs_m, float delta,
    arma::mat Calpha, arma::mat alpha, arma::mat A_m)
{
    int n_tau = tau_m_old_other_dim.n_rows;
    int m_tau = tau_m_old_other_dim.n_cols;
    arma::mat tau_new(n_tau, m_tau);

    if (d == 1)
    {
        //     tau_new <-t(matrix(.xlogy(self$Cpi[[1]][, m],
        //                                self$pi [[m]] [[d]], eps = NULL),
        //                         self$Q[d], self$n[[1]][m])) +
        //                   ((self$nonNAs [[m]]) * self$A [[m]]) % * %
        //                       t(self$Cpi[[2]][, m] * t(self$tau [[m]][[2]])) % * %
        //                       t(.logit(self$Calpha * self$delta[m] * self$alpha,
        //                                eps = 1e-9)) +
        //                   (self$nonNAs [[m]]) % * %
        //                       t(self$Cpi[[2]][, m] * t(self$tau [[m]][[2]])) % * %
        //                       t(.log(1 - self$Calpha * self$alpha * self$delta[m],
        //                              eps = 1e-9))
        tau_new = (Cpi_1_m % log(pi_m_d)).t() +
                  (nonNAs_m % A_m) * (Cpi_2_m % tau_m_old_other_dim.t()).t() *
                      (logit(delta * Calpha % alpha)).t() +
                  nonNAs_m * (Cpi_2_m % tau_m_old_other_dim.t()).t() *
                      (log(1 - Calpha % alpha * delta));

        return tau_new;
    }
    else
    {
        // Assuming d == 2
        return tau_new;
    }
}