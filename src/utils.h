// A collection of common variables and functions

// Copyright (c) 2023-2024 
// Code by Hongbo Yao (hongbo.yao@outlook.com) and Zhengyong Ren 
// (renzhengyong@csu.edu.cn), with contributions by Zijun Zuo 
// (zzjdqqh@163.com) and Chaojian Chen (chaojian.chen@lmu.de)
// https://github.com/hongbo-yao

#ifndef _UTILS_H
#define _UTILS_H

#include <complex>

typedef std::complex<double>  Dcomplex;
static const Dcomplex   II    = Dcomplex(0.,1.);
static const double     pi    = 3.1415926535897932384626433832795;
static const double     mu0   = 4.0*pi*1e-7;
static const double     R0    = 6371.2;  // Earth's mean radius in km


#endif // _UTILS_H