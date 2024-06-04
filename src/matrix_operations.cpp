#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
SEXP matvecmult_eigen(const Eigen::Map<Eigen::MatrixXd> A,
                      const Eigen::Map<Eigen::VectorXd> b,
                      int n_cores){

  Eigen::setNbThreads(n_cores);
  Eigen::VectorXd C = A * b;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP matTvecmult_eigen(const Eigen::Map<Eigen::MatrixXd> A,
                       const Eigen::Map<Eigen::VectorXd> b,
                       int n_cores){

  Eigen::setNbThreads(n_cores);
  Eigen::VectorXd C = A.transpose() * b;
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
SEXP matmatTmult_eigen(const Eigen::Map<Eigen::MatrixXd> A,
                       const Eigen::Map<Eigen::MatrixXd> B,
                       int n_cores){

  Eigen::setNbThreads(n_cores);
  Eigen::MatrixXd C = A*B.transpose();
  return Rcpp::wrap(C);
}


