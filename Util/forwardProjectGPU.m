function projection = forwardProjectGPU(psf,Xguess)
%% Forward projection from 3D volume to 2D projetion (using CPU)
%% Input:
% @psf: the phase space PSF  
% @Xguess: the projection obtained from "forward project"
%% Output:
% projection: estimate volume  
%
%   [1]  JIAMIN WU, ZHI LU and DONG JIANG etc,
%        Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale
%        Cell, 2021. 
%
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

projection=gpuArray.zeros(size(Xguess,1),size(Xguess,2),'single');

[ra ca]=size(Xguess(:,:,1));
[rb cb]=size(psf(:,:,1));
r = ra+rb-1;c=ca+cb-1; p1 = (r-ra)/2;
a1 = gpuArray.zeros(r,r,'single');
b1 = gpuArray.zeros(r,r,'single');
for z=1:size(Xguess,3)
  %%  
    a1(1:ra,1:ca) = Xguess(:,:,z) ;
    b1(1:rb,1:cb) = psf(:,:,z) ;
    clear con1;
    con1 = ifft2(fft2(a1).*fft2(b1));
    projection = projection + real(con1(p1+1:r-p1,p1+1:r-p1));
end
end
