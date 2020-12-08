% 3D reconstruciton with DAO (for single-frame)
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
addpath('./Solver/');
addpath('./Util/');

% Preparameters
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
maxIter=20; %% the maximum iteration number 
simPSF=0; %% using simulated PSF (on/off)
GPUcompute=1; %% GPU accelerator (on/off)
AO=1;     %% digital adaptive optics (on/off)
if simPSF==1
    load('../PSFmatrix/psf_matrix.mat'); %% the filepath of PSF
else
    load('../PSFmatrix/psf_M40_NA1.0_zmin-15u_zmax_30u.mat','psf'); %% the filepath of PSF
end
weight=squeeze(sum(sum(sum(psf,1),2),5))./sum(psf(:));
weight=weight-min(weight(:));
weight=weight./max(weight(:)).*0.2;
for u=1:Nnum
    for v=1:Nnum
        if (u-round(Nnum/2))^2+(v-round(Nnum/2))^2>16 
            weight(u,v)=0;
        end
    end
end
load('../Data/04_Timeweighted/exampleData/20191211Testis.mat'); %% the filepath of timeweighted multiplexed phase-space data

% Initialization
WDF=imresize(WDF,[size(WDF,1)*Nnum/Nshift,size(WDF,2)*Nnum/Nshift]);
Xguess=ones(size(WDF,1),size(WDF,2),size(psf,5));
Xguess=Xguess./sum(Xguess(:)).*sum(WDF(:))./(size(WDF,3)*size(WDF,4));


% 3D deconvolution with DAO
Xguess = deconvRL(maxIter, Xguess,WDF, psf,weight,AO,GPUcompute);

% save high-resolution reconstructed volume
imwriteTFSK(gather(Xguess),'../Data/05_Recon/exampleData/Recon3D.tif');
