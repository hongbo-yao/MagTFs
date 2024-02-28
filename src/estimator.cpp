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

#include "estimator.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>
#include <cassert>
#include <cfloat> // DBL_MAX
#include <omp.h>

Estimator::Estimator(Param& param_): param(&param_)
{

}



Estimator::~Estimator()
{

}



void Estimator::estimate_TF_from_time_series()
{
   // The number of input and output channels
   int n_output_channels =  param->output_data_list[0][0].size();
   int n_input_channels =  param->input_data_list[0][0].size();

   // Initialization
   std::vector<std::vector<std::vector<int>>> estimated_flag(param->n_periods);
   TF.resize(param->n_periods);
   delta_TF.resize(param->n_periods);
   coh2.resize(param->n_periods);
   coh2_mult.resize(param->n_periods);
   for (int f=0; f<param->n_periods; f++)
   {
      estimated_flag[f].resize(n_output_channels);
      TF[f].resize(n_output_channels);
      delta_TF[f].resize(n_output_channels);
      coh2[f].resize(n_output_channels);
      coh2_mult[f].resize(n_output_channels);
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      {
         estimated_flag[f][output_channel_id].resize(n_input_channels,0);
         TF[f][output_channel_id].resize(n_input_channels);
         delta_TF[f][output_channel_id].resize(n_input_channels);
         coh2[f][output_channel_id].resize(n_input_channels);
      }
   }

   // Loop all periods
   for (int f=0; f<param->n_periods; f++)
   {
      /* Compute the number of segments for this period */
      double T = param->periods[f];
      // section length KT
      double KT = param->K*T;
      // number of data points of this segment
      int n_data = std::round(KT/param->sample_interval);
      std::vector<double> hamming = hamming_window(n_data);

      // length of time series (all pieces)
      // NOTE: all pieces are independent!
      std::vector<double> L_list(param->n_real_input_channel_pieces);
      std::vector<int> N_segments_max_list(param->n_real_input_channel_pieces);
      std::vector<int> N_segments_max_list_idx(param->n_real_input_channel_pieces,0);
      // number of maximum segments
      int N_segments_max = 0;
      for (int i=0; i<param->n_real_input_channel_pieces; i++)
      {
         L_list[i] = param->Nt_list[i]*param->sample_interval;
         N_segments_max_list[i] = std::floor((L_list[i]-2*KT)/((1-param->overlap)*KT))+3;
         N_segments_max += N_segments_max_list[i];
         N_segments_max_list_idx[i+1] = N_segments_max;
      }
      if (N_segments_max < 4*n_input_channels) break;

      /* Loop all pieces and segments to assemble the least square linear system */
      // frequency domain input signals
      MatrixXcd A_max(N_segments_max,n_input_channels); A_max.setZero();
      // frequency domain output signals
      std::vector<VectorXcd> b_max(n_output_channels);
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      {
         b_max[output_channel_id].resize(N_segments_max);
         b_max[output_channel_id].setZero();
      }
      std::vector<int> segment_has_gap(N_segments_max,0);
      int n_gap_segments = 0;
      // Loop all pieces
      for (int piece_id=0; piece_id<param->n_real_input_channel_pieces; piece_id++)
      {
         std::vector<std::vector<Dcomplex>> input_data = param->input_data_list[piece_id];
         std::vector<std::vector<Dcomplex>> output_data = param->output_data_list[piece_id];
         // Loop all segments
         #pragma omp parallel for reduction(+:n_gap_segments)
         for (int s=0; s<N_segments_max_list[piece_id]; s++)
         {
            int linear_system_idx = N_segments_max_list_idx[piece_id]+s;
            // std::cout << linear_system_idx << "\n";
            int idx = s*(int)(KT*(1-param->overlap)/param->sample_interval);
            bool delete_segment = false;

            /* Assemble b using output channel data */
            for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
            {
               // step 1: prepare the segments
               std::vector<Dcomplex> output_seg(n_data);
               std::vector<double> output_gap_idx;
               std::vector<Dcomplex> output_seg_no_gap;
               std::vector<double> output_seg_no_gap_idx;
               for (int t=0; t<n_data; t++)
               {
                  output_seg[t] = output_data[idx+t][output_channel_id];
                  if (output_seg[t]!=param->gap)
                  {
                     output_seg_no_gap.push_back(output_data[idx+t][output_channel_id]);
                     output_seg_no_gap_idx.push_back(t);
                  }
                  else output_gap_idx.push_back(t);
               }
               // Interpolating the gap data
               int n_output_gaps = n_data-output_seg_no_gap_idx.size();
               double output_gap_percent = (double)n_output_gaps/(double)n_data;
               if (n_output_gaps!=0)
               {
                  if (output_gap_percent<=param->gap_interp_floor)
                  {
                     std::vector<Dcomplex> output_seg_gap_value = 
                        interp1(output_seg_no_gap_idx, output_seg_no_gap, output_gap_idx);
                     for (int i=0; i<output_seg_gap_value.size(); i++)
                     {
                        output_seg[(int)output_gap_idx[i]] = output_seg_gap_value[i];
                     }
                  }
                  // To much gaps in this segment, delete it
                  else 
                  {
                     delete_segment = true;
                     segment_has_gap[linear_system_idx] = 1;
                     if (output_channel_id==0) n_gap_segments++;
                  }
               }
               if (delete_segment) continue;

               // step 2: detrend
               detrend(output_seg);

               // step 3: apply a window function to reduce sidelobes
               for (int t=0; t<n_data; t++)
               {
                  output_seg[t] *= hamming[t];
               }

               // step 4: fourier tramsform
               for (int t=0; t<n_data; t++)
               {
                  Dcomplex ft = std::exp(-II*2.0*pi/T*(double)t*param->sample_interval);
                  b_max[output_channel_id][linear_system_idx] += output_seg[t]*ft;
               }
            }
            if (delete_segment) continue;

            /* Assemble A using input channel data */
            for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
            {
               // step 1: prepare the segments
               std::vector<Dcomplex> input_seg(n_data);
               std::vector<double> input_gap_idx; // using double for interp1
               std::vector<Dcomplex> input_seg_no_gap;
               std::vector<double> input_seg_no_gap_idx;
               for (int t=0; t<n_data; t++)
               {
                  input_seg[t] = input_data[idx+t][input_channel_id];
                  if (input_seg[t]!=param->gap)
                  {
                     input_seg_no_gap.push_back(input_seg[t]);
                     input_seg_no_gap_idx.push_back(t);
                  }
                  else input_gap_idx.push_back(t);
               }
               // Interpolating the gap data
               int n_input_gaps = n_data-input_seg_no_gap_idx.size();
               double input_gap_percent = (double)n_input_gaps/(double)n_data;
               if (n_input_gaps!=0)
               {
                  // Interpolating the gap data
                  if (input_gap_percent<=param->gap_interp_floor)
                  {
                     std::vector<Dcomplex> input_seg_gap_value = 
                        interp1(input_seg_no_gap_idx, input_seg_no_gap, input_gap_idx);
                     for (int i=0; i<input_seg_gap_value.size(); i++)
                     {
                        input_seg[(int)input_gap_idx[i]] = input_seg_gap_value[i];
                     }
                  }
                  // To much gaps in this segment, delete it
                  else 
                  {
                     delete_segment = true;
                     segment_has_gap[linear_system_idx] = 1;
                     n_gap_segments++;
                  }
               }
               if (delete_segment) break;

               // step 2: detrend
               detrend(input_seg);

               // step 3: apply a window function to reduce sidelobes
               for (int t=0; t<n_data; t++)
               {
                  input_seg[t] *= hamming[t];
               }

               // step 4: fourier tramsform
               for (int t=0; t<n_data; t++)
               {
                  Dcomplex ft = std::exp(-II*2.0*pi/T*(double)t*param->sample_interval);
                  A_max(linear_system_idx,input_channel_id) += input_seg[t]*ft;
               }
            }
         }
      }

      // Delete segments with gap
      int N_segments = N_segments_max-n_gap_segments;
      if (N_segments < 4*n_input_channels) break;

      // least-square linear system Ax=b, with Eigen for dense linear operation
      // frequency domain input signals and TFs
      MatrixXcd A(N_segments,n_input_channels); A.setZero();
      // frequency domain output signals, TFs, and squared coherence
      std::vector<VectorXcd> b(n_output_channels);
      std::vector<VectorXcd> x(n_output_channels);
      std::vector<VectorXd> coh2_(n_output_channels);
      std::vector<double> coh2_mult_(n_output_channels,0.);
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      {
         b[output_channel_id].resize(N_segments);
         b[output_channel_id].setZero();
         x[output_channel_id].resize(N_segments);
         x[output_channel_id].setZero();
         coh2_[output_channel_id].resize(n_input_channels);
         coh2_[output_channel_id].setZero();
      }
      int idx = 0;
      for (int i=0; i<N_segments_max; i++)
      {
         if (!segment_has_gap[i])
         {
            for (int j=0; j<n_input_channels; j++) A(idx,j) = A_max(i,j);
            for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
            {
               b[output_channel_id][idx] = b_max[output_channel_id][i];
            }
            idx++;
         }
      }
      A_max.resize(0,0);
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      {
         b_max[output_channel_id].resize(0);
      }
      b_max.resize(0);

      // Print infomation
      if (T>86400)
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T/86400 << " days, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }
      else if (T>3600)
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T/3600 << std::setw(10) << " hours, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }
      else
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T << std::setw(10) << " seconds, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }

      // Loop all output channels
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      {
         // Solve the least square system 
         irls(A, b[output_channel_id], x[output_channel_id], coh2_[output_channel_id], coh2_mult_[output_channel_id]);
         coh2_mult[f][output_channel_id] = coh2_mult_[output_channel_id];
         for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
         {
            TF[f][output_channel_id][input_channel_id] = x[output_channel_id][input_channel_id];
            coh2[f][output_channel_id][input_channel_id] = coh2_[output_channel_id][input_channel_id];
            estimated_flag[f][output_channel_id][input_channel_id] = 1;
         }

         // Compute the data uncertainty using the jackknife method (Chave and Thomson, 1989)
         // This part was referred to Püthe, C. (2015), pages 29-30 and 
         // Semenov and Kuvshinov (2012), page 969
         std::vector<std::vector<Dcomplex>> TF_leave_one_out(N_segments);
         int rows = N_segments-1;
         int cols = n_input_channels;
         double TF_average_re[n_input_channels] = {0.}, TF_average_im[n_input_channels] = {0.};
         #pragma omp parallel for reduction(+:TF_average_re,TF_average_im)
         for (int s=0; s<N_segments; s++)
         {
            MatrixXcd A_temp = A;
            VectorXcd b_temp = b[output_channel_id];
            VectorXcd x_temp(N_segments-1); x_temp.setZero();
            VectorXd coh2_temp(n_input_channels);
            double coh2_mult_temp;

            // delete the s-th row, faster than loop all rows and cols
            A_temp.block(s,0,rows-s,cols) = A_temp.block(s+1,0,rows-s,cols);
            A_temp.conservativeResize(rows,cols);
            b_temp.segment(s,rows-s) = b_temp.segment(s+1,rows-s);
            b_temp.conservativeResize(rows);

            // solve
            irls(A_temp, b_temp, x_temp, coh2_temp, coh2_mult_temp);
            TF_leave_one_out[s].resize(n_input_channels);
            for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
            {
               TF_leave_one_out[s][input_channel_id] = x_temp[input_channel_id];
               TF_average_re[input_channel_id] += std::real(TF_leave_one_out[s][input_channel_id]);
               TF_average_im[input_channel_id] += std::imag(TF_leave_one_out[s][input_channel_id]);
            }
         }

         for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
         {
            Dcomplex TF_average(TF_average_re[input_channel_id],TF_average_im[input_channel_id]);
            TF_average /= N_segments;

            // jackknife estimation of transfer function
            Dcomplex TF_jk =  (double)N_segments*TF[f][output_channel_id][input_channel_id]-((double)N_segments-1)*TF_average;

            // jackknife estimation of variance
            double C_jk = 0;
            #pragma omp parallel for reduction(+:C_jk)
            for (int i=0; i<N_segments; i++)
            {
               C_jk += std::pow(std::abs(TF_leave_one_out[i][input_channel_id]-TF_average),2);
            }
            C_jk *= (double)(N_segments-n_input_channels)/(double)N_segments;

            // data uncertainty
            // delta_TF[f][output_channel_id][input_channel_id] = std::sqrt(C_jk)*tinv(1.0-(1.0-param->confidence_level)/2.0,N_segments-n_input_channels);
            delta_TF[f][output_channel_id][input_channel_id] = std::sqrt(C_jk);
         }
      }
   }

   // // Check uncertainty, if uncertainty is greater than both real 
   // // and imaginary parts, then delete this period
   // for (int i=0; i<param->n_periods; i++)
   // {
   //    for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
   //    { 
   //       for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
   //       {
   //          if (estimated_flag[i][output_channel_id][input_channel_id])
   //          {
   //             if (delta_TF[i][output_channel_id][input_channel_id]>std::abs(TF[i][output_channel_id][input_channel_id].real()) && delta_TF[i][output_channel_id][input_channel_id]>std::abs(TF[i][output_channel_id][input_channel_id].imag()))
   //             {
   //                estimated_flag[i][output_channel_id][input_channel_id] = 0;
   //             }
   //          }
   //       }
   //    }
   // }

   // Output TFs, data uncertainty, and coh2
   std::string TF_file = param->output_file_prefix + ".TF";
   std::ofstream out_stream(TF_file);
   out_stream << "# period_id, output_channel_id, input_channel_id, period, TF_re, TF_im, TF_std_err, coh2, coh2_mult\n";
   for (int i=0; i<param->n_periods; i++)
   {
      for (int output_channel_id=0; output_channel_id<n_output_channels; output_channel_id++)
      { 
         for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
         {
            if (estimated_flag[i][output_channel_id][input_channel_id])
            {
               out_stream  << std::setw(6) << std::setprecision(6) << i+1
                           << std::setw(6) << std::setprecision(6) << output_channel_id+1
                           << std::setw(6) << std::setprecision(6) << input_channel_id+1
                           << std::setw(16) << std::setprecision(8) << param->periods[i]
                           << std::setw(16) << std::setprecision(6) << TF[i][output_channel_id][input_channel_id].real()
                           << std::setw(16) << std::setprecision(6) << TF[i][output_channel_id][input_channel_id].imag()
                           << std::setw(16) << std::setprecision(6) << delta_TF[i][output_channel_id][input_channel_id] 
                           << std::setw(16) << std::setprecision(6) << coh2[i][output_channel_id][input_channel_id]
                           << std::setw(16) << std::setprecision(6) << coh2_mult[i][output_channel_id] <<"\n";
            }
         }
      }
   }
}



void Estimator::estimate_TF_from_field_spectra()
{
   int n_periods = param->field_spectra_list[0].period_id.size();
   std::vector<std::vector<int>> estimated_flag(n_periods);
   for (int f=0; f<n_periods; f++)
   {
      int n_input_channels = param->source_list[0][f].SH_degree_list.size();
      estimated_flag[f].resize(n_input_channels,0);
   }

   std::vector<std::vector<Dcomplex>> spectra_TF(n_periods);
   std::vector<std::vector<double>> spectra_delta_TF(n_periods);
   std::vector<std::vector<double>> spectra_coh2(n_periods);
   std::vector<double> spectra_coh2_mult(n_periods);
   for (int f=0; f<n_periods; f++)
   {
      int n_input_channels = param->source_list[0][f].SH_degree_list.size();
      double T = param->source_list[0][f].period;
      int N_segments = param->field_spectra_list.size();
      if (N_segments < 4*n_input_channels) break;
      // std::cout << n_input_channels << "\n";
      spectra_TF[f].resize(n_input_channels);
      spectra_delta_TF[f].resize(n_input_channels);
      spectra_coh2[f].resize(n_input_channels);

      // Print infomation
      if (T>86400)
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T/86400 << " days, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }
      else if (T>3600)
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T/3600 << std::setw(10) << " hours, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }
      else
      {
         std::cout << "Estimating TFs at period #" << std::setw(2) << f+1 
                   << ": " << std::setw(8) << std::setiosflags(std::ios::fixed) 
                   << std::setprecision(2) << T << std::setw(10) << " seconds, " 
                   << std::setw(7) << N_segments << " segments" << "\n";
      }

      // Least-square linear system Ax=b, with Eigen for dense linear operation
      // frequency domain input signals
      MatrixXcd A(N_segments,n_input_channels); A.setZero();
      // frequency domain output signals
      VectorXcd b(N_segments); b.setZero();
      // frequency domain TFs
      VectorXcd x(N_segments); x.setZero();
      // squared coherence
      VectorXd coh2_(n_input_channels);
      double coh2_mult_;
      for (int i=0; i<N_segments; i++)
      {
         b[i] = param->field_spectra_list[i].field_spectra[f];
         SH_source source = param->source_list[i][f];
         for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
         {
            A(i,input_channel_id) = source.e_nm_list[input_channel_id];
         }
      }

      // Solve the least square system 
      irls(A, b, x, coh2_,coh2_mult_);
      spectra_coh2_mult[f] = coh2_mult_;
      for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
      {
         spectra_TF[f][input_channel_id] = x[input_channel_id];
         spectra_coh2[f][input_channel_id] = coh2_[input_channel_id];
         estimated_flag[f][input_channel_id] = 1;
      }

      // Compute the data uncertainty using the jackknife method (Chave and Thomson, 1989)
      // This part was referred to Püthe, C. (2015), pages 29-30 and 
      // Semenov and Kuvshinov (2012), page 969
      std::vector<std::vector<Dcomplex>> TF_leave_one_out(N_segments);
      int rows = N_segments-1;
      int cols = n_input_channels;
      double TF_average_re[n_input_channels] = {0.}, TF_average_im[n_input_channels] = {0.};
      #pragma omp parallel for reduction(+:TF_average_re,TF_average_im)
      for (int s=0; s<N_segments; s++)
      {
         MatrixXcd A_temp = A;
         VectorXcd b_temp = b;
         VectorXcd x_temp(N_segments-1); x_temp.setZero();
         VectorXd coh2_temp(n_input_channels);
         double coh2_mult_temp;

         // delete the s-th row, faster than loop all rows and cols
         A_temp.block(s,0,rows-s,cols) = A_temp.block(s+1,0,rows-s,cols);
         A_temp.conservativeResize(rows,cols);
         b_temp.segment(s,rows-s) = b_temp.segment(s+1,rows-s);
         b_temp.conservativeResize(rows);

         // solve
         irls(A_temp, b_temp, x_temp, coh2_temp, coh2_mult_temp);
         TF_leave_one_out[s].resize(n_input_channels);
         for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
         {
            TF_leave_one_out[s][input_channel_id] = x_temp[input_channel_id];
            TF_average_re[input_channel_id] += std::real(TF_leave_one_out[s][input_channel_id]);
            TF_average_im[input_channel_id] += std::imag(TF_leave_one_out[s][input_channel_id]);
         }
      }

      for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
      {
         Dcomplex TF_average(TF_average_re[input_channel_id],TF_average_im[input_channel_id]);
         TF_average /= N_segments;

         // jackknife estimation of transfer function
         Dcomplex TF_jk =  (double)N_segments*spectra_TF[f][input_channel_id]-((double)N_segments-1)*TF_average;

         // jackknife estimation of variance
         double C_jk = 0;
         #pragma omp parallel for reduction(+:C_jk)
         for (int i=0; i<N_segments; i++)
         {
            C_jk += std::pow(std::abs(TF_leave_one_out[i][input_channel_id]-TF_average),2);
         }
         C_jk *= (double)(N_segments-n_input_channels)/(double)N_segments;

         // data uncertainty
         // spectra_delta_TF[f][input_channel_id] = std::sqrt(C_jk)*tinv(1.0-(1.0-param->confidence_level)/2.0,N_segments-n_input_channels);
         spectra_delta_TF[f][input_channel_id] = std::sqrt(C_jk);
      }
   }

   // Output TFs, data uncertainty, and spectra_coh2
   std::stringstream TFfile;
   TFfile << param->output_file_prefix << ".TF";
   std::ofstream TF_stream(TFfile.str().c_str());
   TF_stream << "# period_id, period, n, m, Tnm_re, Tnm_im, uncertainty, coh2, coh2_mult\n";
   for (int f=0; f<n_periods; f++)
   {
      int period_id = param->source_list[0][f].period_id;
      double T = param->source_list[0][f].period;
      SH_source source = param->source_list[0][f];
      int n_input_channels = source.SH_degree_list.size();
      for (int input_channel_id=0; input_channel_id<n_input_channels; input_channel_id++)
      {
         if (estimated_flag[f][input_channel_id])
         {
            TF_stream   << std::setw(4) << period_id 
                        << std::setw(8) << std::setprecision(8) << T
                        << std::setw(4) << source.SH_degree_list[input_channel_id]
                        << std::setw(4) << source.SH_order_list[input_channel_id]
                        << std::setw(14) << std::setprecision(6) << spectra_TF[f][input_channel_id].real()
                        << std::setw(14) << std::setprecision(6) << spectra_TF[f][input_channel_id].imag() 
                        << std::setw(14) << std::setprecision(6) << spectra_delta_TF[f][input_channel_id]
                        << std::setw(14) << std::setprecision(6) << spectra_coh2[f][input_channel_id]
                        << std::setw(14) << std::setprecision(6) << spectra_coh2_mult[f] << "\n";
         }
      }
   }
   TF_stream.close();
}



int Estimator::find_nearest_neighbor_idx(std::vector<double> x, double value)
{
   int idx = 0;
   double dist = DBL_MAX;
   for (int i=0; i<x.size(); i++)
   {
      double dist2 = value-x[i];
      if (dist2>=0 && dist2<dist)
      {
         idx = i;
         dist = dist2;
      }
   }
   return idx;
}



std::vector<Dcomplex> Estimator::interp1(std::vector<double> x, 
                                         std::vector<Dcomplex> y, 
                                         std::vector<double> x_new)
{
   std::vector<Dcomplex> y_new(x_new.size());
   std::vector<Dcomplex> slope(x.size());
   std::vector<Dcomplex> intercept(x.size());
   // interpolation
   for (int i=0; i<x.size()-1; i++)
   {
      slope[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
      intercept[i] = y[i]-slope[i]*x[i];
   }
   // extrapolation
   slope[x.size()-1] = slope[x.size()-2];
   intercept[x.size()-1] = intercept[x.size()-2];

   for (int i=0; i<x_new.size(); i++)
   {
      int idx = find_nearest_neighbor_idx(x, x_new[i]);
      y_new[i] = slope[idx]*x_new[i]+intercept[idx];
   }
   return y_new;
}



std::vector<double> Estimator::hamming_window(int N)
{
   std::vector<double> hamming(N);
   double a0 = 0.54;
   for (int n=0; n<N; n++)
   {
      hamming[n] = a0-(1-a0)*std::cos(2*pi*n/(N-1));
   }
   return hamming;
}



void Estimator::detrend_real_data(std::vector<double>& y)
{
   int m = y.size();
   std::vector<double> x(m);
   for (int i=0; i<m; i++)
   {
      x[i] = i;
   }

   double xmean = 0;
   double ymean = 0;
   for (int i=0; i<m; i++)
   {
      xmean += x[i];
      ymean += y[i];
   }
   xmean /= (double)m;
   ymean /= (double)m;

   // calculate Covariance
   double temp = 0;
   for (int i=0; i<m; i++) 
   {
      temp += x[i]*y[i];
   }
   double Sxy = temp/(double)m-xmean*ymean;
   temp = 0;
   for (int i=0; i<m; i++) 
   {
      temp += x[i]*x[i];
   }
   double Sxx = temp/(double)m-xmean*xmean;

   double grad = Sxy/Sxx;
   double yint = -grad*xmean+ymean;

   // remove linear trend
   for (int i=0; i<m; i++)
   {
      y[i] = y[i]-(grad*(double)i+yint);
   }
}


void Estimator::detrend(std::vector<Dcomplex>& y)
{
   int m = y.size();
   std::vector<double> y_real(m);
   std::vector<double> y_imag(m);
   for (int i=0; i<m; i++)
   {
      y_real[i] = y[i].real();
      y_imag[i] = y[i].imag();
   }
   detrend_real_data(y_real);
   detrend_real_data(y_imag);
   for (int i=0; i<m; i++)
   {
      y[i] = Dcomplex(y_real[i],y_imag[i]);
   }
}



void Estimator::irls(MatrixXcd G, VectorXcd d, VectorXcd& m, 
                     VectorXd& coh2_temp, double& coh2_mult_temp, 
                     int maxit, double rtol)
{
   const int Nd = G.rows();
   const int n_input_channels = G.cols();
   // complex conjugate transpose, do not use transpose()!
   MatrixXcd Gt = G.adjoint(); 
   RowVectorXcd dt = d.adjoint();
   MatrixXcd GtG = Gt*G; // pxp
   VectorXcd Gtb = Gt*d; // px1
   Dcomplex dtd = dt*d; // 1x1
   VectorXd e;
   RowVectorXd et;
   double ete;

   // Least square
   if (param->least_square_method=="ls")
   {
      m = GtG.colPivHouseholderQr().solve(Gtb);
      e = (d-G*m).cwiseAbs();
      ete = e.transpose()*e;
   }
   // IRLS with Huber weight
   else
   {
      // the first iteration, weight W=1
      m = GtG.colPivHouseholderQr().solve(Gtb);
      VectorXcd m_old = m;
      e = (d-G*m).cwiseAbs();
      ete = e.transpose()*e;
      double r = std::sqrt(ete/Nd);
      VectorXd W(Nd);
      for (int i=1; i<maxit; i++)
      {
         double Wsum = 0.;
         for (int j=0; j<Nd; j++)
         {
            double temp = 1.5*r/e[j];
            W[j] = (temp < 1.0) ? temp : 1.0;
            Wsum += W[j];
         }
         MatrixXd Wmat = W.replicate(1,n_input_channels);
         VectorXcd temp_vec = W.array()*d.array();
         GtG = Gt*(Wmat.cwiseProduct(G)); // Gt*W*G
         Gtb = Gt*temp_vec; // Gt*W*d
         dtd = dt*temp_vec; // dt*W*d
         m = GtG.colPivHouseholderQr().solve(Gtb);

         // Compute standard derivation
         e = (d-G*m).cwiseAbs(); // |Ax-d|
         VectorXd We = W.array()*e.array();
         ete = e.transpose()*We;
         r = std::sqrt(ete/Wsum);

         // Check convergence
         double tol = (m-m_old).squaredNorm()/(m_old.squaredNorm());
         m_old = m;
         if (tol < rtol) break;
      }
   }

   // Compute squared coherence coh2, eq (6) in Semenov and Kuvshinov (2012)
   for (int i=0; i<n_input_channels; i++)
   {
      coh2_temp[i] = std::real(std::pow(std::abs(Gtb(i)),2)/(GtG(i,i)*dtd));
   }
   coh2_mult_temp = std::real(1.0-ete/dtd);
}