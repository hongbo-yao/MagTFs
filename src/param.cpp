// Geomagnetic transfer function estimation parameters

// Copyright (c) 2023-2024 
// Code by Hongbo Yao (hongbo.yao@outlook.com) and Zhengyong Ren 
// (renzhengyong@csu.edu.cn), with contributions by Zijun Zuo 
// (zzjdqqh@163.com) and Chaojian Chen (chaojian.chen@lmu.de)
// https://github.com/hongbo-yao

#include "param.h"
#include <fstream>
#include <cassert>
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>

Param::Param(char* config_file)
{
   read_input_parameters(config_file);
   // using time series as input
   if (segment_dir_list_file.empty()) 
   {
      read_time_series();
      if (period_type=="auto") set_periods();
   }
   // using field spectra as input
   else 
   {
      read_field_spectra();
   }
}



Param::~Param()
{

}



void Param::read_input_parameters(char* config_file)
{
   std::ifstream in_stream(config_file);
   assert(in_stream.good());
   int n_real_input_channel_readin = 0;
   int n_imag_input_channel_readin = 0;
   int n_real_output_channel_readin = 0;
   int n_imag_output_channel_readin = 0;
   int n_periods_readin = 0;
   std::string line;
   std::string key,value,comment;
   while (std::getline(in_stream, line))
   {
      // skip line that starts with '#' or empty
      if (*(line.begin()) == '#' || line == "") continue;
      int loc = line.find(':');
      if(loc != std::string::npos)
      {
         // obtain key
         key = line.substr(0,loc);
         // obtain value
         int comment_loc = line.find('#');
         if(comment_loc == std::string::npos) // no comments
            value = line.substr(loc+1,line.size());
         else
            value = line.substr(loc+1,comment_loc);

         // set value
         if (key.find("Number of real input channel pieces")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_real_input_channel_pieces;
         }
         if (key.find("Number of imag input channel pieces")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_imag_input_channel_pieces;
         }
         if (key.find("Number of real output channel pieces")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_real_output_channel_pieces;
         }
         if (key.find("Number of imag output channel pieces")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_imag_output_channel_pieces;
         }
         else if (key.find("Output file prefix")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> output_file_prefix;
         }
         else if (key.find("Gap id")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> gap_id;
         }
         else if (key.find("Gap interp floor")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> gap_interp_floor;
         }
         else if (key.find("Sample interval")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> sample_interval;
         }
         else if (key.find("Section length multiple")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> K;
         }
         else if (key.find("Overlap")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> overlap;
         }
         else if (key.find("Least square method")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> least_square_method;
         }
         else if (key.find("Period type")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> period_type;
         }
         else if (key.find("Number of periods")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> n_periods;
         }
         else if (key.find("Minimum period")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> period_min;
         }
         else if (key.find("Maximum period")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> period_max;
         }
         else if (key.find("Segment directory")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> segment_dir;
         }
         else if (key.find("segment_dir_list_file")!=std::string::npos)
         {
            std::stringstream value_sstr(value);
            value_sstr >> segment_dir_list_file;
         }
      }

      if (loc==std::string::npos && n_real_input_channel_readin<n_real_input_channel_pieces)
      {
         std::stringstream sstr(line);
         std::string data_file;
         sstr >> data_file;
         real_input_channel_file_list.push_back(data_file);
         n_real_input_channel_readin++;
      }

      if (loc==std::string::npos && n_imag_input_channel_readin<n_imag_input_channel_pieces)
      {
         std::stringstream sstr(line);
         std::string data_file;
         sstr >> data_file;
         imag_input_channel_file_list.push_back(data_file);
         n_imag_input_channel_readin++;
      }

      if (loc==std::string::npos && n_real_output_channel_readin<n_real_output_channel_pieces)
      {
         std::stringstream sstr(line);
         std::string data_file;
         sstr >> data_file;
         real_output_channel_file_list.push_back(data_file);
         n_real_output_channel_readin++;
      }

      if (loc==std::string::npos && n_imag_output_channel_readin<n_imag_output_channel_pieces)
      {
         std::stringstream sstr(line);
         std::string data_file;
         sstr >> data_file;
         imag_output_channel_file_list.push_back(data_file);
         n_imag_output_channel_readin++;
      }

      if (loc==std::string::npos && period_type=="given" && n_periods_readin<n_periods)
      {
         std::stringstream sstr(line);
         double temp_period;
         sstr >> temp_period; 
         periods.push_back(temp_period);
         n_periods_readin++;
      }
   }

   gap = Dcomplex(gap_id,0.);
   assert(n_real_input_channel_pieces==n_real_output_channel_pieces);
   if (n_imag_input_channel_pieces>0)
   {
      assert(n_imag_input_channel_pieces==n_real_input_channel_pieces);
      assert(n_imag_input_channel_pieces==n_imag_output_channel_pieces);
      gap = Dcomplex(gap_id,gap_id);
   }

   station_name = output_file_prefix;
   if (!segment_dir_list_file.empty() && station_name=="observed")
   {
      std::cout << "Wrong output_file_prefix! Please set it as station name for global-to-local TFs!\n";
      std::abort();
   }
}



std::vector<std::vector<double>> Param::load(std::string data_file)
{
   std::ifstream in_stream(data_file);
   assert(in_stream.good());
   std::vector<std::vector<double>> data;
   std::string line;
   while (std::getline(in_stream, line))
   {
      if (*(line.begin()) == '#' || line == "") continue;
      std::vector<double> row;
      std::stringstream sstr(line);
      double tmp;
      while(sstr>>tmp)
      {
         row.push_back(tmp);
      }
      data.push_back(row);
   }
   in_stream.close();
   return data;
}



void Param::read_time_series()
{
   input_data_list.resize(n_real_input_channel_pieces);
   output_data_list.resize(n_real_output_channel_pieces);
   Nt_list.resize(n_real_output_channel_pieces);
   for (int piece=0; piece<n_real_input_channel_pieces; piece++)
   {
      // input data
      std::string real_input_channel_file = real_input_channel_file_list[piece];   
      std::vector<std::vector<double>> real_input_data = load(real_input_channel_file);
      std::vector<std::vector<Dcomplex>> input_data(real_input_data.size());
      for (int i=0; i<real_input_data.size(); i++)
      {
         input_data[i].resize(real_input_data[i].size());
         for (int j=0; j<real_input_data[i].size(); j++)
         {
            input_data[i][j] = Dcomplex(real_input_data[i][j],0.);
         }
      }
      if (n_imag_input_channel_pieces>0)
      {
         std::string imag_input_channel_file = imag_input_channel_file_list[piece];
         std::vector<std::vector<double>> imag_input_data = load(imag_input_channel_file);
         for (int i=0; i<imag_input_data.size(); i++)
         {
            for (int j=0; j<imag_input_data[i].size(); j++)
            {
               if (real_input_data[i][j]==gap_id) input_data[i][j] = Dcomplex(gap_id,gap_id);
               else input_data[i][j] += Dcomplex(0.,imag_input_data[i][j]);
            }
         }
      }
      input_data_list[piece] = input_data;
      // for (int i=0; i<input_data.size(); i++)
      // {
      //    for (int j=0; j<input_data[i].size(); j++)
      //    {
      //       std::cout << input_data[i][j] << "\t";
      //    }
      //    std::cout << "\n";
      // }

      // output data
      std::string real_output_channel_file = real_output_channel_file_list[piece];   
      std::vector<std::vector<double>> real_output_data = load(real_output_channel_file);
      std::vector<std::vector<Dcomplex>> output_data(real_output_data.size());
      for (int i=0; i<real_output_data.size(); i++)
      {
         output_data[i].resize(real_output_data[i].size());
         for (int j=0; j<real_output_data[i].size(); j++)
         {
            output_data[i][j] = Dcomplex(real_output_data[i][j],0.);
         }
      }
      if (n_imag_output_channel_pieces>0)
      {
         std::string imag_output_channel_file = imag_output_channel_file_list[piece];
         std::vector<std::vector<double>> imag_output_data = load(imag_output_channel_file);
         for (int i=0; i<imag_output_data.size(); i++)
         {
            for (int j=0; j<imag_output_data[i].size(); j++)
            {
               if (real_output_data[i][j]==gap_id) output_data[i][j] = Dcomplex(gap_id,gap_id);
               else output_data[i][j] += Dcomplex(0.,imag_output_data[i][j]);
            }
         }
      }
      output_data_list[piece] = output_data;
      // for (int i=0; i<output_data.size(); i++)
      // {
      //    for (int j=0; j<output_data[i].size(); j++)
      //    {
      //       std::cout << output_data[i][j] << "\t";
      //    }
      //    std::cout << "\n";
      // }

      int Nt_in = input_data_list[piece].size();
      int Nt_out = output_data_list[piece].size();
      assert(Nt_in==Nt_out);
      Nt_list[piece] = Nt_in;
   }
}



void Param::set_periods()
{
   periods.resize(n_periods);
   double a = std::log10(period_min);
   double b = std::log10(period_max);
   double step = (b-a)/(n_periods-1);
   for (int i=0; i<n_periods; i++)
   {
      double period = a+i*step;
		periods[i] = std::ceil(std::pow(10, period));
      // std::cout << periods[i] << "\n";
   }
}



void Param::read_field_spectra()
{
   std::vector<std::string> segment_dir_list;
   std::ifstream segment_dir_list_in_stream(segment_dir_list_file);
   assert(segment_dir_list_in_stream.good());
   std::string line;
   while(std::getline(segment_dir_list_in_stream,line))
   {
      if (*(line.begin()) == '#' || line == "") continue;
      segment_dir_list.push_back(line);
   }
   segment_dir_list_in_stream.close();
   // for (int i=0; i<segment_dir_list.size(); i++)
   // {
   //    std::cout << segment_dir_list[i] << "\n";
   // }

   for (int i=0; i<segment_dir_list.size(); i++)
   {
      // read field spectra
      Field_spectra tmp_field_spectra;
      std::ostringstream field_spectra_filename;
      field_spectra_filename << segment_dir << "/" << segment_dir_list[i] << "/" 
                             << station_name << segment_dir_list[i] << ".Zspectra";
      // std::cout << field_spectra_filename.str() << "\n";
      std::ifstream field_spectra_in_stream(field_spectra_filename.str().c_str());
      assert(field_spectra_in_stream.good());
      while(std::getline(field_spectra_in_stream,line))
      {
         if (*(line.begin()) == '#' || line == "") continue;
         std::stringstream sstr(line);
         int tmp_period_id;
         double tmp_period, Z_re, Z_im;
         sstr >> tmp_period_id >> tmp_period >> Z_re >> Z_im;
         tmp_field_spectra.period_id.push_back(tmp_period_id);
         tmp_field_spectra.period.push_back(tmp_period);
         tmp_field_spectra.field_spectra.push_back(Dcomplex(Z_re,Z_im));
      }
      field_spectra_in_stream.close();
      field_spectra_list.push_back(tmp_field_spectra);

      // read source
      std::vector<SH_source> tmp_source;
      for (int j=0; j<tmp_field_spectra.period_id.size(); j++)
      {
         std::ostringstream filename;
         filename << segment_dir << "/" << segment_dir_list[i] << "/period" 
                  << tmp_field_spectra.period_id[j] << "_source.txt";
         // std::cout << filename.str() << "\n";
         std::ifstream in_stream(filename.str().c_str());
         SH_source source;
         source.period_id = tmp_field_spectra.period_id[j];
         source.period = tmp_field_spectra.period[j];
         while(std::getline(in_stream,line))
         {
            if (*(line.begin()) == '#' || line == "") continue;
            std::stringstream sstr(line);
            double n, m, e_nm_re, e_nm_im;
            sstr >> n >> m >> e_nm_re >> e_nm_im;
            source.SH_degree_list.push_back(n);
            source.SH_order_list.push_back(m);
            source.e_nm_list.push_back(Dcomplex(e_nm_re, e_nm_im));
         }
         in_stream.close();
         tmp_source.push_back(source);
      }
      source_list.push_back(tmp_source);
   }
}