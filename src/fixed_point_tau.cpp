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
    return log(x / (1.0 - x));
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

// [[Rcpp::export]]
arma::mat fixed_point_tau(
    int d, arma::mat tau_m_old_other_dim,
    arma::mat Cpi_1_m, arma::mat Cpi_2_m,
    arma::mat pi_m_d, int Q_d, int n_d_m,
    arma::mat nonNAs_m, float delta,
    arma::mat Calpha, arma::mat alpha, arma::mat A_m)
{
    int n_tau = tau_m_old_other_dim.n_rows;
    int m_tau = tau_m_old_other_dim.n_cols;
    arma::mat tau_new(n_tau, m_tau);

    if (d == 1)
    {
        tau_new = log(Cpi_1_m % pi_m_d) + (nonNAs_m % A_m) * (Cpi_2_m % tau_m_old_other_dim) * (logit(delta * alpha)).t() + nonNAs_m * (Cpi_2_m % tau_m_old_other_dim) * (log(ones(size(alpha)) - delta * Calpha % alpha)).t();
        return tau_new;
    }
    else if (d == 2)
    {
        tau_new = log(Cpi_2_m % pi_m_d) + (nonNAs_m % A_m).t() * (Cpi_1_m % tau_m_old_other_dim) * logit(delta * alpha) + nonNAs_m.t() * (Cpi_1_m % tau_m_old_other_dim) * log(ones(size(alpha)) - delta * Calpha % alpha);
        return tau_new;
    }
    else
    {
        return tau_new;
    }
}