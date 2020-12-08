function Backprojection = backwardProjectACC(psf,projection)
%% Backprojection from 2D projetion to 3D volume (using CPU)
%% Input:
% @psf: the phase space PSF  
% @projection: the projection obtained from "forward project". 
%% Output:
% Backprojection: estimate volume  
%
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020
psf=rot90(psf,2);
Backprojection=zeros(size(projection,1),size(projection,2),size(psf,3));
for z=1:size(psf,3)
    Backprojection(:,:,z)=conv2(projection,psf(:,:,z),'same');
end

