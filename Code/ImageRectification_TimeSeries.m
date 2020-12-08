% Image rectification of scanning light field data (for time-series)
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
filepath='../Data/01_Raw/exampleData_Timeseries/20191107Neutrophil.tif'; %% the filepath of raw scanning light field data
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
Nx=30; %% half of microlens number in the first dimension
Ny=30; %% half of microlens number in the second dimension
centerX=1088; %% central coordinate in the first dimension
centerY=724; %% central coordinate in the second dimension

% scanning order
if Nshift==3
    index1=[3,2,1,1,1,2,3,3,2];
    index2=[1,1,1,2,3,3,3,2,2];
end

% rectification for time series
for frame=0:99
    sLF=zeros( (2*Nx+1)*Nnum,(2*Ny+1)*Nnum,Nshift,Nshift );
    for i=1:Nshift*Nshift
        sLF_slice=double(imread(filepath, i+frame));
        ii = mod((i+frame)-1,Nshift*Nshift)+1;
        sLF(:,:,index2(ii),index1(ii))=sLF_slice(centerX-Nnum*Nx-fix(Nnum/2):centerX+Nnum*Nx+fix(Nnum/2),centerY-Nnum*Ny-fix(Nnum/2):centerY+Nnum*Ny+fix(Nnum/2));
    end
    % save rectificated scanning light field data
    save(['../Data/02_Rectify/exampleData_Timeseries/20191107Neutrophil_frame',num2str(frame),'.mat'],'sLF','-v7.3');
end
