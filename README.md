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

1. Unpack the package
2. Include subdirectory in your Matlab path
3. Run main_for_static.m or main_for_timelapse.m for static or time-lapse data.
   a). The scanning Lightfield Data referred to Fig.1F (the imaging of neutrophils in mouse liver), Fig.S1F (the simulation of fluorescent beads) and Fig.S5C (the imaging of the testis slice) in our paper are saved in "Data" which can be used for test. 
   b). MAT file of the experimental calibrated point spread function is located in dir "PSF". Readers can also generate the similar ideal PSFs by running "main_computePSF.m" (based on wave optics theory[1]), which is located in dir "PSFcalculation".
* Download the test scanning light filed images (neutrophils in mouse liver, the testis slice) and experimental PSFs from the following link.
https://drive.google.com/drive/folders/101IHbAApPF-Z734UtjDOZHEZwtBleQgC?usp=sharing

----------------
Main modules description
----------------
1. main_for_static_testis.m: Pre-processes and 3D reconstruciton with DAO (for single-frame static data) (>=70 GB memory)
2. main_for_static_beads.m: Pre-processes and 3D reconstruciton with DAO (for single-frame static data) (>=16 GB memory)
3. main_for_timelapse_miceliver.m: Pre-processes and 3D reconstruciton with DAO (for time-lapse data) (>=40 GB memory)
* Pre-processes including ImageRectification, Realignment, Timeweighted algorithm and 3D deconvolution with DAO, are involved in main_for_static.m and main_for_timelapse.m.

----------------
IMPORTANT NOTE 
----------------
Should you have any questions regarding this code and the corresponding results, please contact Zhi Lu (luz18@mails.tsinghua.edu.cn).

Reference:
1.  R. Prevedel, Y.-G. Yoon, M. Hoffmann, N. Pak, G. Wetzstein, S. Kato, T. Schr?del, R. Raskar, M. Zimmer, E. S. Boyden, and A. Vaziri, 
     "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy," Nature Methods, 2014, 11(7): 727-730.
2.  Lu Z, Wu J, Qiao H, et al. "Phase-space deconvolution for light field microscopy," Optics express, 2019, 27(13): 18131-18145.
