# MagTF
MagTF is a C++ code for estimating geomagnetic transfer functions from observatory and satellite magnetic data. Supported transfer functions include:

- Geomagnetic depth sounding global scalar Q-responses (using time series of external and internal SH coefficients as input)

- Geomagnetic depth sounding global matrix Q-responses (using time series of external and internal SH coefficients as input)

- Geomagnetic depth sounding local C-responses (using time series of geomagnetic observatory data as input)

- Geomagnetic depth sounding global-to-local transfer functions (using time series or field spectra of external SH coefficients and geomagnetic observatory data as input)

- Magnetotelluric scalar impedance and tipper (using time series of geomagnetic observatory data as input)

Main features include:

- Parallel computing with C++ and OpenMP

- Can handle time series with gaps

- Can handle time series with multiple pieces

- Can deal with user-defined discrete periods

- Flexible input file format and easy-to-use

### 1. Prerequisite
Eigen is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms (https://eigen.tuxfamily.org/index.php?title=Main_Page).

MagTF uses Eigen for some linear algebra operations. After downloading the Eigen library, please put it into contrib/ directory. A copy of Eigen also can be found in contrib/ directory

### 2. Usage
The command to use MagTF is 

```
MagTF config_file
```

### 3. Citation
MagTF was developed and used in the following papers:

#### Geomagnetic observatory data and local C-responses
- Hongbo Yao, Zhengyong Ren*, Jingtian Tang, Keke Zhang. A multi-resolution finite-element approach for global electromagnetic induction modeling with application to southeast China coastal geomagnetic observatory studies. Journal of Geophysical Research: Solid Earth, 2022, 127(8), e2022JB024659. https://doi.org/10.1029/2022JB024659

- Hongbo Yao, Zhengyong Ren*, Jingtian Tang, Rongwen Guo, Jiayong Yan. Trans-dimensional Bayesian joint inversion of magnetotelluric and geomagnetic depth sounding responses to constrain mantle electrical discontinuities. Geophysical Journal International, 2023, https://doi.org/10.1093/gji/ggad029

#### Satellite magnetic data and global Q-responses
- Hongbo Yao, Zhengyong Ren*, Kejia Pan, Jingtian Tang, Keke Zhang. A global mantle conductivity model derived from 8 years of Swarm satellite magnetic data. Earth and Planetary Physics, 2023, 7(1), 49-56. https://doi.org/10.26464/epp2023011

#### Geomagnetic observatory data and global-to-local transfer functions
- Please follow my researchgate profile https://www.researchgate.net/profile/Hongbo-Yao-2

#### Satellite magnetic data and matrix Q-responses
- Please follow my researchgate profile https://www.researchgate.net/profile/Hongbo-Yao-2

If using this code, a citation to our papers is appreciated.

### 4. References
The theory behind MagTF code is based on the following references:

- Semenov, A., & Kuvshinov, A. (2012). Global 3-D imaging of mantle conductivity based on inversion of 
observatory C-responses-II. Data analysis and results. Geophysical Journal International, 191(3), 965–992. 
https://doi.org/10.1111/j.1365-246X.2012.05665.x

- Püthe, C. (2015). Interpretation of global EM induction data from ground, sea and space. New response functions, 
inversion schemes and conductivity models [ETH Zurich]. https://doi.org/10.3929/ethz-a-010531597

- Chave, A. D., & Thomson, D. J. (1989). Some comments on magnetotelluric response function estimation. 
Journal of Geophysical Research, 94(B10). https://doi.org/10.1029/jb094ib10p14215

### 5. Contact
- Hongbo Yao, Email: hongbo.yao@outlook.com or yaohongbogeo@foxmail.com

- Zhengyong Ren, Email: renzhengyong@csu.edu.cn
