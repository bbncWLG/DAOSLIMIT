% Timeweighted algorithm (for single-frame): adjust the weight in different
%                    ...scanning sampling points to remove motion artefacts
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
addpath('./Util/');

% Preparameters
Nshift=3; %% the sampling points of a single scanning period
Nnum=13; %% number of virtual pixels
time_weight_index=1; %% timeweighted coefficient, ranging from 0 to 1

% scanning order
if Nshift==3
    a1=[1,1,1,2,3,3,3,2,2];
    a2=[3,2,1,1,1,2,3,3,2];
end

% Data processing
load('../Data/03_Realign/exampleData/20191211Testis.mat'); %% the filepath of multiplexed phase-space
index=zeros( size(WDF)  );
now1=mod(0,Nshift*Nshift); %% the scanning position at the moment
index1=[a1(now1+1:end),a1(1:now1)]; %% the scanning order at the moment
index2=[a2(now1+1:end),a2(1:now1)]; %% the scanning order at the moment
WDF=double(WDF);

% weight distribuiton according to time and scanning positions
for time_aver=1:1:Nshift*Nshift+1  
    if time_aver == Nshift*Nshift+1
        time_mask_wigner=ones(size(WDF,1),size(WDF,2),size(WDF,3),size(WDF,4));
    else
        time_mask_wigner=zeros(size(WDF,1),size(WDF,2),size(WDF,3),size(WDF,4));
        time_mask_wigner(Nshift+1-index1(time_aver):Nshift:end,Nshift+1-index2(time_aver):Nshift:end,:,: )=1;
    end
    WDF_t=WDF.*time_mask_wigner;
    if time_aver ~= Nshift*Nshift+1
        index(WDF_t ~= 0) =time_weight_index^(abs(time_aver-round(Nshift^2/2)));
    end
end

% TimeWeighted
WDF_TW =zeros(size(WDF,1),size(WDF,2),Nnum,Nnum);
for i = 1:Nnum
    for j = 1:Nnum
        datanew = Timeweighted_interp(WDF_t(:,:,i,j),index(:,:,i,j),a1,a2,Nshift);
        datanew(isnan(datanew))=0;
        WDF_TW(:,:,i,j)=gather(datanew);
    end
end
WDF_TW(WDF_TW<0)=0;
WDF=WDF_TW(Nshift+1:end-Nshift,Nshift+1:end-Nshift,:,:);

% save timeweighted multiplexed phase-space data
save('../Data/04_Timeweighted/exampleData/20191211Testis.mat','WDF','-v7.3');


