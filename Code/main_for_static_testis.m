% Pre-processes and 3D reconstruciton with DAO (for single-frame static data)
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
filepath='../Data/01_Raw/exampleData_Static/20191211Testis.tif'; %% the filepath of raw scanning light field data
Nx=50; %% half of used microlens number in the first dimension
Ny=50; %% half of uesd microlens number in the second dimension
centerX=1064; %% central coordinate in the first dimension
centerY=1123; %% central coordinate in the second dimension
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
maxIter=10; %% the maximum iteration number 
ExperimentalPSF=1; %% using experimental PSF or ideal PSF
GPUcompute=1; %% GPU accelerator (on/off)
DAO=1;     %% digital adaptive optics (on/off)
time_weight_index=1; %% timeweighted coefficient, ranging from 0 to 1
frame=0; %% frame=0 if the data is static

% scanning order
if Nshift==3
    index1=[3,2,1,1,1,2,3,3,2];
    index2=[1,1,1,2,3,3,3,2,2];
end
if ExperimentalPSF==0
    load('../PSFmatrix/Ideal_psf_M40_NA1.0_zmin-15u_zmax30u.mat'); %% the filepath of ideal PSF
else
    load('../PSFmatrix/Experimental_psf_M40_NA1.0_zmin-15u_zmax30u.mat','psf'); %% the filepath of experimental PSF
end
weight=squeeze(sum(sum(sum(psf,1),2),5))./sum(psf(:));
weight=weight-min(weight(:));
weight=weight./max(weight(:)).*0.8;
for u=1:Nnum
    for v=1:Nnum
        if (u-round(Nnum/2))^2+(v-round(Nnum/2))^2>16 
            weight(u,v)=0;
        end
    end
end

% rectification
sLF=zeros( (2*Nx+1)*Nnum,(2*Ny+1)*Nnum,Nshift,Nshift );
for i=1:Nshift*Nshift        
    sLF_slice=double(imread(filepath, i));
    sLF(:,:,index2(i),index1(i))=sLF_slice(centerX-Nnum*Nx-fix(Nnum/2):centerX+Nnum*Nx+fix(Nnum/2),centerY-Nnum*Ny-fix(Nnum/2):centerY+Nnum*Ny+fix(Nnum/2));
end

% Pixel realignment
multiWDF=zeros(Nnum,Nnum,size(sLF,1)/Nnum,size(sLF,2)/Nnum,Nshift,Nshift); %% multiplexed phase-space
for i=1:Nnum
    for j=1:Nnum
        for a=1:size(sLF,1)/Nnum
            for b=1:size(sLF,2)/Nnum
                multiWDF(i,j,a,b,:,:)=squeeze(  sLF(  (a-1)*Nnum+i,(b-1)*Nnum+j,:,:  )  );
            end
        end
    end
end
WDF=zeros(  size(sLF,1)/Nnum*Nshift,size(sLF,2)/Nnum*Nshift,Nnum,Nnum  ); % multiplexed phase-space
for a=1:size(sLF,1)/Nnum
    for c=1:Nshift
        x=Nshift*a+1-c;
        for b=1:size(sLF,2)/Nnum
            for d=1:Nshift
                y=Nshift*b+1-d;
                WDF(x,y,:,:)=squeeze(multiWDF(:,:,a,b,c,d));
            end
        end
    end
end


% Time-weighted
WDF=time_weighted(WDF,time_weight_index,index1,index2,Nshift,Nnum,frame);

% Initialization
WDF=imresize(WDF,[size(WDF,1)*Nnum/Nshift,size(WDF,2)*Nnum/Nshift]);
Xguess=ones(size(WDF,1),size(WDF,2),size(psf,5));
Xguess=Xguess./sum(Xguess(:)).*sum(WDF(:))./(size(WDF,3)*size(WDF,4));

% 3D deconvolution with DAO
Xguess = deconvRL(maxIter, Xguess,WDF, psf,weight,DAO,GPUcompute);

% save high-resolution reconstructed volume
imwriteTFSK(single(gather(Xguess(21:end-20,21:end-20,11:end-10))),['../Data/02_Recon/exampleData_Static/Recon3D_20191211Testis.tif']);
