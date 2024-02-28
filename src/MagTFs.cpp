// Estimating geomagnetic transfer functions using the section-averaing 
// approach and Iterative Reweighted Least Squares (IRLS) with Huber weight.
// Data uncertainty is estimated using the jackknife method.

// Copyright (c) 2023-2024 
// Code by Hongbo Yao (hongbo.yao@outlook.com) and Zhengyong Ren 
// (renzhengyong@csu.edu.cn), with contributions by Zijun Zuo 
// (zzjdqqh@163.com) and Chaojian Chen (chaojian.chen@lmu.de)
// https://github.com/hongbo-yao

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

#include <iostream>
#include <ctime>
#include <omp.h>
#include "param.h"
#include "estimator.h"

int main(int argc, char* argv[])
{
   if (argc < 2)
   {
      std::cout << "Usage: " << argv[0] << " config_filename\n"; 
      return 1;    
   }
   double tic = omp_get_wtime();
   Param param(argv[1]);
   Estimator estimator(param);
   if (param.segment_dir_list_file.empty()) estimator.estimate_TF_from_time_series();
   else estimator.estimate_TF_from_field_spectra();
   double toc = omp_get_wtime();
   std::cout << "Run time: " << toc-tic << " (s)\n";
   return 0;
}