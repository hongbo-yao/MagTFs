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

The details of MagTFs is formally described in the following manuscript:

> Zhengyong Ren, Zijun Zuo, Hongbo Yao*, Chaojian Chen*, Linan Xu, Jingtian Tang, Keke Zhang. MagTFs: A tool for estimating multiple magnetic transfer functions to constrain Earth's electrical conductivity structure. Computers & Geosciences, 2025, 195, 105769. https://doi.org/10.1016/j.cageo.2024.105769

Some functions of this tool were developed and used in these studies:

> Hongbo Yao, Zhengyong Ren*, Jingtian Tang, Keke Zhang. A multi-resolution finite-element approach for global electromagnetic induction modeling with application to southeast China coastal geomagnetic observatory studies. Journal of Geophysical Research: Solid Earth, 2022, 127(8), e2022JB024659. https://doi.org/10.1029/2022JB024659

> Hongbo Yao, Zhengyong Ren*, Jingtian Tang, Rongwen Guo, Jiayong Yan. Trans-dimensional Bayesian joint inversion of magnetotelluric and geomagnetic depth sounding responses to constrain mantle electrical discontinuities. Geophysical Journal International, 2023, 233(3), 1821-1846. https://doi.org/10.1093/gji/ggad029

> Hongbo Yao, Zhengyong Ren*, Kejia Pan, Jingtian Tang, Keke Zhang. A global mantle conductivity model derived from 8 years of Swarm satellite magnetic data. Earth and Planetary Physics, 2023, 7(1), 49-56. https://doi.org/10.26464/epp2024047

> Zhengyong Ren, Yifei Xie, Chaojian Chen*, Hongbo Yao*, Jingtian Tang, Keke Zhang. New insights into Earth's mantle conductivity and water distribution from the Macau Science Satellite-1 data. Earth and Planetary Physics, 2025, accepted.

If this helps you, a citation to our papers is appreciated.

### 5. Authors
Code by Hongbo Yao (hongbo.yao@outlook.com), with contributions from Zhengyong Ren (renzhengyong@csu.edu.cn), Zijun Zuo (zzjdqqh@163.com) and Chaojian Chen (chenchaojian@csu.edu.cn)