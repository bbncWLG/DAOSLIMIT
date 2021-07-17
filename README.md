# DAOSLIMIT
Title:      Matlab code for "Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale"
Version:    1.0 
Copyright:  2021, JIAMIN WU, ZHI LU, DONG JIANG, JINGTAO FAN, LI YU, QIONGHAI DAI

Edited based on the reference [1][2].


Matlab code for "Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale"
==========================================================

This package contains the implementations of the algorithms, including pixel realignment, mutual iterative tomography with DAO, time-weighted algorithm, and time-loop reconstructions for time-lapse data, which are described in detail in the paper: 
JIAMIN WU, ZHI LU and DONG JIANG etc, "Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale". Cell 2021.

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) in an academic publication.

For algorithmic details, please refer to our paper.

Additional example data and experimental PSFs can be downloaded with the following link:
https://drive.google.com/drive/folders/101IHbAApPF-Z734UtjDOZHEZwtBleQgC?usp=sharing

----------------
How to use
----------------
The code is tested in MATLAB 2018b (64bit) under the MS Windows 10 64bit version with an Intel i9-9980XE CPU, NVIDIA GeForce RTX 2080 Ti GPU and 128GB RAM.

1. Unpack the package
2. Include subdirectory in your Matlab path
3. Run the .m files with the prefix of "main" for static or time-lapse data.

   a). The raw data of the scanning light field should be placed in the folder of "Data". We have provided several example data for test, including, time-lapse video of neutrophils in mouse liver (63x/1.4NA oil immersion objective), simulated fluorescent beads (63x/1.4NA oil immersion objective), the autofluorescence of a fixed testis slice (40x/1.0NA water immersion objective) and the membrane-labelled zebrafish embryos (63x/1.4NA oil immersion objective). Some of the data are too large for GitHub. So we upload them in the google drive described before. 
   
   b). The PSF files in the format of .mat should be placed in the folder of "PSF". Some of the PSFs are too large, especially for experimental PSFs. So we have uploaded them in the same link described before. Readers can also generate the simulated ideal PSFs by running "main_computePSF.m" (based on wave optics theory and phase-space theory [1][2]), which is located in the folder of "PSFcalculation".

* Download more test images and experimental PSFs from the following link.
https://drive.google.com/drive/folders/101IHbAApPF-Z734UtjDOZHEZwtBleQgC?usp=sharing


----------------
Main modules description
----------------
1. main_for_static_testis.m: Pre-processes and 3D reconstruciton with DAO (for single-frame static data) (>=70 GB memory)
2. main_for_static_beads.m: Pre-processes and 3D reconstruciton with DAO (for single-frame static data) (>=16 GB memory)
3. main_for_static_zebraembryos.m: Pre-processes and 3D reconstruciton with DAO (for single-frame static data) (>=32 GB memory)
3. main_for_timelapse_miceliver.m: Pre-processes and 3D reconstruciton with DAO, time-weighted algorithm and time-loop algorithm (for time-lapse data) (>=40 GB memory)
* Pre-processes including ImageRectification and pixel realignment are involved in main_for_static.m and main_for_timelapse.m.

----------------
IMPORTANT NOTE 
---------------- 
Should you have any questions regarding this code and the corresponding results, please contact Zhi Lu (luz18@mails.tsinghua.edu.cn).

Reference:
1.  R. Prevedel, Y.-G. Yoon, M. Hoffmann, N. Pak, G. Wetzstein, S. Kato, T. Schr?del, R. Raskar, M. Zimmer, E. S. Boyden, and A. Vaziri, 
     "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy," Nature Methods, 2014, 11(7): 727-730.
2.  Lu Z, Wu J, Qiao H, et al. "Phase-space deconvolution for light field microscopy," Optics express, 2019, 27(13): 18131-18145.
