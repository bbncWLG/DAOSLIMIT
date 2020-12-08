% Pixel realignment (for time-series): from scanning light field data to multiplexed phase-space
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
Nnum=13; %% number of virtual pixels
Nshift=3; %% the sampling points of a single scanning period
for frame=0:99 %% time-series
load(['../Data/02_Rectify/exampleData_Timeseries/20191107Neutrophil_frame',num2str(frame),'.mat']); %% the filepath of rectificated scanning light field data

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
WDF=zeros(  size(sLF,1)/Nnum*Nshift,size(sLF,2)/Nnum*Nshift,Nnum,Nnum  ); %% multiplexed phase-space
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

% save multiplexed phase-space
save(['../Data/03_Realign/exampleData_TimeSeries/20191107Neutrophil_frame',num2str(frame),'.mat'],'WDF','-v7.3');
end