function projection = forwardProjectACC(psf,Xguess)
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
%    Contact: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

projection=zeros(size(Xguess,1),size(Xguess,2));
for z=1:size(Xguess,3)
    projection=projection+conv2(Xguess(:,:,z),psf(:,:,z),'same');
end

end

