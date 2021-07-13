function Backprojection = backwardProjectGPU(psf,projection)
%% Backprojection from 2D projetion to 3D volume (using GPU)
%% Input:
% @psf: the phase space PSF  
% @projection: the projection obtained from "forward project". 
%% Output:
% Backprojection: estimate volume  
%
%   [1]  JIAMIN WU, ZHI LU and DONG JIANG etc,
%        Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale
%        Cell, 2021. 
%
%    Contact: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

psf=rot90(psf,2);
Backprojection=gpuArray.zeros(size(projection,1),size(projection,2),size(psf,3),'single');
[ra ca]=size(projection);
[rb cb]=size(psf(:,:,1));
r = ra+rb-1;c=ca+cb-1; p1 = (r-ra)/2;
b1 = gpuArray.zeros(r,r,'single');
a1 = gpuArray.zeros(r,r,'single');
for z=1:size(psf,3)
    a1(1:ra,1:ca) = projection(:,:) ;
    b1(1:rb,1:cb) = psf(:,:,z) ;
    clear con1;
    con1 = ifft2(fft2(a1) .* fft2(b1));
    Backprojection(:,:,z) = real(con1(p1+1:r-p1,p1+1:r-p1));
    
end
end
