function xout = zeroPad(xin, zeroImage)
%% Robert's code
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.
xout = zeroImage;
xout( 1:size(xin,1), 1:size(xin,2) ) = xin;