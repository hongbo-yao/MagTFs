// Geomagnetic transfer function estimation parameters

// Copyright (c) 2023-2024 
// Code by Hongbo Yao (hongbo.yao@outlook.com) and Zhengyong Ren 
// (renzhengyong@csu.edu.cn), with contributions by Zijun Zuo 
// (zzjdqqh@163.com) and Chaojian Chen (chaojian.chen@lmu.de)
// https://github.com/hongbo-yao

#ifndef _PARAM_H
#define _PARAM_H

#include <string>
#include <vector>
#include "utils.h"

// Spherical harmonic sources for one period
struct SH_source
{
   int period_id;
   double period;
   std::vector<double> SH_degree_list; // n
   std::vector<double> SH_order_list;  // m
   std::vector<Dcomplex> e_nm_list;
};

// Field spectra for one segment directory
struct Field_spectra
{
   std::vector<int> period_id;
   std::vector<double> period;
   std::vector<Dcomplex> field_spectra;
};


class Param
{
public:
   // Input and output channels with multiple pieces
   // The imaginary parts can be elimated if the input and output channels are real-valued
   int n_real_input_channel_pieces = 0;
   int n_imag_input_channel_pieces = 0;
   int n_real_output_channel_pieces = 0;
   int n_imag_output_channel_pieces = 0;
   std::vector<std::string> real_input_channel_file_list;
   std::vector<std::string> imag_input_channel_file_list;
   std::vector<std::string> real_output_channel_file_list;
   std::vector<std::string> imag_output_channel_file_list;
   // n_data_pieces x Nt x n_input_channels
   std::vector<std::vector<std::vector<Dcomplex>>> input_data_list;
   // n_data_pieces x Nt x n_output_channels
   std::vector<std::vector<std::vector<Dcomplex>>> output_data_list;
   // Number of time points for each data piece
   std::vector<int> Nt_list;

   // Prefix of output filename
   std::string output_file_prefix = "example";

   // Gap marker in the data
   double gap_id;
   Dcomplex gap;

   // If gap percent<gap_interp_floor, we apply linear interpolation
   double gap_interp_floor = 0.02;

   // Sample interval in seconds
   double sample_interval;

   // KT is the section length for period T
   int K = 3;

   // Overlap percentage of adjacent sections
   double overlap = 0.5;

   // Least square method used to solve the least square system
   // 'ls' denotes regular least square method
   // 'irls' denotes iterative reweighted least squares
   std::string least_square_method = "ls";

   // Periods
   std::string period_type = "";
   int n_periods;
   double period_min;
   double period_max;
   std::vector<double> periods;


   /* Using field spectra as input, mainly designed for estimating Dst and Sq global-to-local transfer functions */
   // In this case, all infomation including periods are included by input data
   // segment_dir_list_file contains all segments directory
   std::string segment_dir_list_file = "";

   // big segment directory, which contains a list of segments directory as provided by segment_dir_list_file
   std::string segment_dir = "";
   
   // Station name, used by read_field_spectra()
   std::string station_name;

   // N_segments x 1
   std::vector<Field_spectra> field_spectra_list;

   // N_segments x n_periods
   std::vector<std::vector<SH_source>> source_list;

public:
   Param(char* config_file);
   ~Param();

public:
   void read_input_parameters(char* config_file);

   // Load a matrix data from file
   std::vector<std::vector<double>> load(std::string data_file);

   // Read data with multiple pieces, woking when using time series as input
   void read_time_series();

   // Read data with field spectra as input
   void read_field_spectra();

   // Set discrete periods according to given min and max periods
   void set_periods();
};

#endif // _PARAM_H