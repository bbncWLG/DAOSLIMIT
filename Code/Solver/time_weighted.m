function WDF=time_weighted(WDF,time_weight_index,index1,index2,Nshift,Nnum,frame)

% Timeweighted algorithm: adjust the weight in different
%                    ...scanning sampling points to remove motion artefacts
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

addpath('./Util/');


% Data processing
WDF=WDF+1e-5;
index=zeros( size(WDF)  );
now1=mod(frame,Nshift*Nshift); %% the scanning position at the moment
now_index1=[index1(now1+1:end),index1(1:now1)]; %% the scanning order at the moment
now_index2=[index2(now1+1:end),index2(1:now1)]; %% the scanning order at the moment

% Weight distribuiton according to time and scanning positions
for time_aver=1:1:Nshift*Nshift+1  
    if time_aver == Nshift*Nshift+1
        time_mask_wigner=ones(size(WDF,1),size(WDF,2),size(WDF,3),size(WDF,4));
    else
        time_mask_wigner=zeros(size(WDF,1),size(WDF,2),size(WDF,3),size(WDF,4));
        time_mask_wigner(Nshift+1-now_index2(time_aver):Nshift:end,Nshift+1-now_index1(time_aver):Nshift:end,:,: )=1;
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
        datanew = Timeweighted_interp(WDF_t(:,:,i,j),index(:,:,i,j),index2,index1,Nshift);
        datanew(isnan(datanew))=0;
        WDF_TW(:,:,i,j)=gather(datanew);
    end
end
WDF_TW(WDF_TW<0)=0;
WDF=WDF_TW(Nshift+1:end-Nshift,Nshift+1:end-Nshift,:,:);



