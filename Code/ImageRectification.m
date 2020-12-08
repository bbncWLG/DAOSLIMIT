% Image rectification of scanning light field data (for single-frame)
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
filepath='../Data/01_Raw/exampleData/20191211Testis.tif'; %% the filepath of raw scanning light field data
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
Nx=50; %% half of microlens number in the first dimension
Ny=50; %% half of microlens number in the second dimension
centerX=998; %% central coordinate in the first dimension
centerY=1043; %% central coordinate in the second dimension

% scanning order
if Nshift==3
    index1=[1,2,3,2,3,1,3,2,1];
    index2=[3,3,3,2,2,2,1,1,1];
end

% rectification
sLF=zeros( (2*Nx+1)*Nnum,(2*Ny+1)*Nnum,Nshift,Nshift );
for i=1:Nshift*Nshift        
    sLF_slice=double(imread(filepath, i));
    sLF(:,:,index2(i),index1(i))=sLF_slice(centerX-Nnum*Nx-fix(Nnum/2):centerX+Nnum*Nx+fix(Nnum/2),centerY-Nnum*Ny-fix(Nnum/2):centerY+Nnum*Ny+fix(Nnum/2));
end

% save rectificated scanning light field data
save('../Data/02_Rectify/exampleData/20191211Testis.mat','sLF','-v7.3');
