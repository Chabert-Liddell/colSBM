#include <RcppArmadillo.h>
using namespace Rcpp;

// // Function to calculate xlogy
// inline NumericMatrix xlogy(const NumericMatrix &x, const NumericMatrix &y, double eps = 1e-9)
// {
//     return log(x + eps) * y;
// }

// Function to calculate log(x+eps)
inline arma::mat log_eps(const arma::mat &x, double eps = 1e-9)
{
    return log(x + eps);
}

// Function to calculate logit
inline arma::mat logit(const arma::mat &x)
{
    return log_eps(x / (1.0 - x));
}

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

// [[Rcpp::depends(RcppArmadillo)]]

// Hoping this allow me to parallelize in R
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::mat fixed_point_tau(
    const int &d, const arma::mat &tau_m_old_other_dim,
    const arma::mat &Cpi_1_m, const arma::mat &Cpi_2_m,
    const arma::mat &pi_m_d,
    const arma::mat &nonNAs_m, const float &delta,
    const arma::mat &Calpha, const arma::mat &alpha, const arma::mat &A_m)
{
    int n_tau = tau_m_old_other_dim.n_rows;
    int m_tau = tau_m_old_other_dim.n_cols;
    arma::mat tau_new(n_tau, m_tau);

    if (d == 1)
    {
        // With the regular log, infinite values are produced on zeros of the 
        // support. Yet they are ignored in the R code and thus non problematic
        tau_new = log_eps(Cpi_1_m % pi_m_d) + (nonNAs_m % A_m) * (Cpi_2_m % tau_m_old_other_dim) * (logit(delta * alpha)).t() + nonNAs_m * (Cpi_2_m % tau_m_old_other_dim) * (log(ones(size(alpha)) - delta * Calpha % alpha)).t();
    }
    else if (d == 2)
    {
        tau_new = log_eps(Cpi_2_m % pi_m_d) + (nonNAs_m % A_m).t() * (Cpi_1_m % tau_m_old_other_dim) * logit(delta * alpha) + nonNAs_m.t() * (Cpi_1_m % tau_m_old_other_dim) * log(ones(size(alpha)) - delta * Calpha % alpha);
        return tau_new;
    }
    else
    {
        tau_new = arma::zeros(size(tau_new));
    }
    if (arma::is_finite(tau_new))
    {
        return tau_new;
    }
    else
    {
        Rcout << "\nError with the tau computation!\n";
        if (d == 1) {
            Rcout << "1st term:\n" << log(Cpi_1_m % pi_m_d);
        } else {
            Rcout << "1st term:\n" << log(Cpi_2_m % pi_m_d);
        }
    }
    // Base return
    return tau_new;
}