
#include "Rcpp.h"
#include "RcppEigen.h"
using namespace Rcpp;
#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "Eigen/LU"
#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include "Eigen/Sparse"

#include<redsvd/redsvd.h>


using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::MatrixXf;
using Eigen::MappedSparseMatrix;
using Eigen::SparseMatrix;


// Note: SEXP AA must be a sparse matrix as in 
// https://cran.r-project.org/web/packages/RcppEigen/vignettes/RcppEigen-Introduction.pdf
// from package Matrix 

// [[Rcpp::export]]
List redSVD(SEXP AA,int num){
	// int num=(int)(INTEGER(nn)[0]);
	
	const MappedSparseMatrix<double> A(as<MappedSparseMatrix<double> >(AA));
	RedSVD::RedSVD<MappedSparseMatrix<double>> svA(A, num);
// return Rcpp::wrap(svA.matrixV());
return List::create(Named("V") = Rcpp::wrap(svA.matrixV()),
Named("U")= Rcpp::wrap(svA.matrixU()),
Named("D")= Rcpp::wrap(svA.singularValues()));
}
