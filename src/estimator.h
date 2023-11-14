// Estimating geomagnetic transfer functions using the section-averaing 
// approach and Iterative Reweighted Least Squares (IRLS) with Huber weight.
// Data uncertainty is estimated using the jackknife method.

// Copyright (c) 2023 Hongbo Yao
// Email: hongbo.yao@outlook.com
// https://github.com/hongbo-yao
// https://www.researchgate.net/profile/Hongbo_Yao2 

// References: 
// Estimating TFs using section-averaging approach:
// [1] Semenov, A., & Kuvshinov, A. (2012). Global 3-D imaging of mantle 
// conductivity based on inversion of observatory C-responses-II. Data 
// analysis and results. Geophysical Journal International, 191(3), 965–992. 
// https://doi.org/10.1111/j.1365-246X.2012.05665.x

// [2] Püthe, C. (2015). Interpretation of global EM induction data from 
// ground, sea and space. New response functions, inversion schemes and 
// conductivity models [ETH Zurich]. https://doi.org/10.3929/ethz-a-010531597

// Estimating data uncertainty using jackknife method:
// [3] Chave, A. D., & Thomson, D. J. (1989). Some comments on magnetotelluric 
// response function estimation. Journal of Geophysical Research, 94(B10). 
// https://doi.org/10.1029/jb094ib10p14215
// see also Püthe, C. (2015).

// IRLS and compute inverse t distribution value:
// [4] Parameter Estimation and Inverse Problems, 3rd edition, 2018 
// by R. Aster, B. Borchers, C. Thurber

#ifndef _ESTIMATOR_H
#define _ESTIMATOR_H

#include "param.h"
#include "utils.h"

#include <Eigen/Dense>  // Linear Algebra Lib.
using namespace Eigen;

class Estimator
{
public:
   Param* param;
   // Nf x n_output_channels x n_output_channels
   std::vector<std::vector<std::vector<Dcomplex>>> TF;
   std::vector<std::vector<std::vector<double>>> delta_TF;
   std::vector<std::vector<std::vector<double>>> coh2;
   std::vector<std::vector<double>> coh2_mult;

public:
   Estimator(Param& param_);
   ~Estimator();

public:
   // Estimate electromagnetic transfer functions using time series as input
   // Cut segments, Apply window, Fourier transform, ...
   void estimate_TF_from_time_series();

   // Estimate electromagnetic transfer functions using field spectra as input
   void estimate_TF_from_field_spectra();

   // Linear interpolation and extrapolation
   int find_nearest_neighbor_idx(std::vector<double> x, double value);
   std::vector<Dcomplex> interp1(std::vector<double> x, 
                                 std::vector<Dcomplex> y, 
                                 std::vector<double> x_new);

   // Compute hamming window values
   std::vector<double> hamming_window(int N);

   // Detrend the linear detrend. This implementation is  
   // modified from the detrend_IP() function of DSPFunctions
   // see https://github.com/lppier/DSPFunctions
   void detrend_real_data(std::vector<double>& y);
   void detrend(std::vector<Dcomplex>& y);

   // Iterative Reweighted Least Squares (IRLS)
   void irls(MatrixXcd G, VectorXcd d, VectorXcd& m, 
             VectorXd& coh2_temp, double& coh2_mult_temp, 
             int maxit=5, double rtol=1e-6);
};

#endif // _ESTIMATOR_H
