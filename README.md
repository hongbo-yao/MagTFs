# MagTFs
MagTFs is a C++ code for estimating multi-source magnetic transfer functions to constrain Earth's electrical conductivity structure. Supported transfer functions include:

- Magnetotelluric tippers

- Sq global-to-local transfer functions

- Dst observatory C-responses

- Dst satellite scalar Q-responses

- Dst satellite matrix Q-responses

Main features include:

- Parallel computing with C++ and OpenMP

- Can handle time series with gaps

- Can handle time series with multiple pieces

- Can deal with user-defined discrete periods

- Flexible input file format and easy-to-use

### 1. Prerequisite
Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms (https://eigen.tuxfamily.org/index.php?title=Main_Page).

MagTFs uses Eigen for some linear algebra operations. After downloading the Eigen library, please put it into contrib/ directory. A copy of Eigen also can be found in contrib/ directory

### 2. Usage
To compile MagTFs, in MagTFs/ directory:

```
make
```

To run MagTFs to estimate different TFs, in examples/03_Dst_observatory_C_responses/FRD/:

```
../../../MagTFs setup.config
```

To show the results, run plot_*.m

### 3. Citation

The background theory of MagTFs is described in the following manuscript:

> Zhengyong Ren, Zijun Zuo, Hongbo Yao*, Chaojian Chen*, Linan Xu, Jingtian Tang, Keke Zhang. MagTFs: A tool for estimating multiple magnetic transfer functions to constrain Earth's electrical conductivity structure. Computers & Geosciences, 2024, accepted.

If this helps you, a citation to our paper is appreciated.

### 5. Authors
Code by Hongbo Yao (hongbo.yao@outlook.com), with contributions from Zhengyong Ren (renzhengyong@csu.edu.cn), Zijun Zuo (zzjdqqh@163.com) and Chaojian Chen (chenchaojian@csu.edu.cn)