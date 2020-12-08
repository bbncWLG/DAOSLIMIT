% 3D reconstruciton with DAO (for single-frame)
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
gpuDevice(2);
addpath('./Solver/');
addpath('./Util/');

% Preparameters
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
maxIter=1; %% the maximum iteration number 
simPSF=0; %% using simulated PSF or captured PSF
GPUcompute=1; %% GPU accelerator (on/off) 
AO=1;     %% digital adaptive optics (on/off)
if simPSF==1 
    load('../PSFmatrix/psf_matrix.mat'); %% the filepath of PSF
else
    load('../PSFmatrix/psf_M63_NA1.4_zmin-12u_zmax_12u.mat','psf'); %% the filepath of PSF
end
weight=squeeze(sum(sum(sum(psf,1),2),5))./sum(psf(:));
weight=weight-min(weight(:));
weight=weight./max(weight(:)).*0.8;

% aperture limit for 
for u=1:Nnum
    for v=1:Nnum
        if (u-round(Nnum/2))^2+(v-round(Nnum/2))^2>16 
            weight(u,v)=0;
        end
    end
end

minFrame=0; %% the started frame
maxFrame=99; %% the end frame
Fstep=1; %% the spacing between adjacent frames

Fcount=minFrame;
for frame = [minFrame: Fstep: maxFrame, maxFrame: -Fstep: minFrame] %% time-series, time-loop
    load(['../Data/04_Timeweighted/exampleData_TimeSeries/20191107Neutrophil_frame',num2str(frame),'.mat']); %% the filepath of timeweighted multiplexed phase-space data
    
    WDF=imresize(WDF,[size(WDF,1)*Nnum/Nshift,size(WDF,2)*Nnum/Nshift]);
    
    % Initialization
    if frame==minFrame && Fcount==minFrame
        Xguess=ones(size(WDF,1),size(WDF,2),size(psf,5));
        Xguess=Xguess./sum(Xguess(:)).*sum(WDF(:))./(size(WDF,3)*size(WDF,4));
    else
        % replaced the uniform initial value with the reconstructed result of the previous frame
        Xguess=0.5 .* (Xguess+ones(size(Xguess)).*sum(Xguess(:))./numel(Xguess)); 
    end
    
    if Fcount==minFrame || Fcount==minFrame+Fstep
        AO = 0; %% DAO off
    else
        AO = 1; %% DAO on
    end
    
    % 3D deconvolution with DAO
    Xguess = deconvRL_TimeSeries(maxIter, Xguess,WDF, psf,weight,AO,GPUcompute);
    
    % save high-resolution reconstructed volume
    if Fcount<=maxFrame
        % would not save the middle results
        imwriteTFSK(gather(Xguess),['../Data/05_Recon/exampleData_TimeSeries/Recon3D_frame',num2str(frame),'.tif']);
    else
        imwriteTFSK(gather(Xguess),['../Data/05_Recon/exampleData_TimeSeries/Timeloop_Recon3D_frame',num2str(frame),'.tif']);
    end
    Fcount=Fcount+Fstep;
end
