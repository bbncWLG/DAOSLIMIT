# DAOSLIMIT
Title:      Matlab code for "3D observation of large-scale subcellular dynamics in vivo at the millisecond scale"  
Author:     JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,  
Version:    1.0 
Copyright:  2019, JIAMIN WU, ZHI LU, DONG JIANG, JINGTAO FAN, LI YU, PENG XI, QIONGHAI DAI

Edited based on the reference [1][2].


Matlab code for "3D observation of large-scale subcellular dynamics in vivo at the millisecond scale"
==========================================================

This package contains an implementation of the DAOSLIMIT algorithm described in the paper: 
JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,  "3D observation of large-scale subcellular dynamics in vivo at the millisecond scale".

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) 
in an academic publication.

For algorithmic details, please refer to our paper.

----------------
How to use
----------------
The code is tested in MATLAB 2018b (64bit) under the MS Windows 10 64bit version with an Intel i9-9980XE CPU, NVIDIA GeForce RTX 2080 Ti GPU and 128GB RAM.

1. unpack the package
2. include code/subdirectory in your Matlab path
3. Try the examples included in this package.
4. Download the required data and PSF through "https://cloud.tsinghua.edu.cn/d/04b591075e184283b606/"
   a). The Lightfield Data referred to Fig.1F (the imaging of neutrophils in mice liver) and Fig.S5C (the imaging of the testis slice) in our paper are saved in "Raw" which can be used for test. 
   b). MAT file of the experimental calibrated point spread function  is located in dir "PSFmatrix". Readers can also generate the similar ideal PSFs by running "computePSF.m".

----------------
Main modules description
----------------
1. computePSF.m: to calculate high-resolution phase-space PSF based on wave optics theory.
2. ImageRectification.m: to rectificate raw scanning LF images for single-frame data
    or ImageRectification_TimeSeries.m: to rectificate raw scanning LF images for time-series data
3. Realign.m: to pixel-realign scanning light field data to multiplexed phase-space for single-frame data
    or Realign_TimeSeries.m: to pixel-realign scanning light field data to multiplexed phase-space for time-series data
4. Timeweighted.m: to adjust the weight in different scanning sampling time-points to remove motion artefacts for single-frame data
    or Timeweighted_TimeSeries.m: to adjust the weight in different scanning sampling time-points to remove motion artefacts for time-series data
5. Reconsruction3D.m: to reconstruct 3D high-resolution volume with DAO for single-frame data
    or Reconsruction3D_TimeSeries.m: to reconstruct 3D high-resolution volume with DAO for time-series data
* Please run the code in the above order.
----------------
IMPORTANT NOTE 
----------------
Should you have any questions regarding this code and the corresponding results, please contact Zhi Lu (luz18@mails.tsinghua.edu.cn).

Reference:
1.  R. Prevedel, Y.-G. Yoon, M. Hoffmann, N. Pak, G. Wetzstein, S. Kato, T. Schr?del, R. Raskar, M. Zimmer, E. S. Boyden, and A. Vaziri, 
     "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy," Nature Methods, 2014, 11(7): 727-730.
2.  Lu Z, Wu J, Qiao H, et al. "Phase-space deconvolution for light field microscopy," Optics express, 2019, 27(13): 18131-18145.
