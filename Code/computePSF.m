
% The Code is created based on the method described in the following paper 
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 
%   [2] Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%       "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy "   
%       in Nature Methods VOL.11 NO.7|July 2014.
% The Code is modified and extended from Robert's code
% 
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

clear;
addpath('./Util/');
%%%%%%%%%%%%%%%%%%%%%%% SIM PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%
NA =        1.4; %% numerical aperture of objective
MLPitch =   100*1e-6; %% pitch of the microlens pitch
Nnum =      13; %% number of virtual pixels
OSR =       3; %% spatial oversampling ratio for computing PSF
n =         1.515; %% refractive index
fml =       2100e-6; %% focal length of the microlens pitch
M =         72.15; %% magnification of the objective lens
lambda =    525*1e-9; %% wavelength
zmax =      10*1e-6; %% the axial location of the highest z-plane with respect to the focal plane
zmin =      -10*1e-6; %% the axial location of the lowest z-plane with respect to the focal plane
zspacing =  0.2*1e-6; %% spacing between adjacent z-planes
eqtol = 1e-10;
k = 2*pi*n/lambda; %% k
k0 = 2*pi*1/lambda; %% k
d = fml;   %% optical distance between the microlens and the sensor
ftl = 165e-3;        %% focal length of tube lens
fobj = ftl/M;  %% focal length of objective lens
fnum_obj = M/(2*NA); %% f-number of objective lens (imaging-side)
fnum_ml = fml/MLPitch; %% f-number of micrl lens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% DEFINE OBJECT SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mod(Nnum,2)==0,
   error(['Nnum should be an odd number']); 
end
pixelPitch = MLPitch/Nnum; %% pitch of virtual pixels

x1objspace = [0]; 
x2objspace = [0];
x3objspace = [zmin:zspacing:zmax]+1e-9; % offset
objspace = ones(length(x1objspace),length(x2objspace),length(x3objspace));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p3max = max(abs(x3objspace));
x1testspace = (pixelPitch/OSR)* [0:1: Nnum*OSR*60];
x2testspace = [0];   
[psfLine] = calcPSFFT(p3max, fobj, NA, x1testspace, pixelPitch/OSR, lambda, d, M, n);
outArea = find(psfLine<0.01);
if isempty(outArea),
   error('Estimated PSF size exceeds the limit');   
end
IMGSIZE_REF = ceil(outArea(1)/(OSR*Nnum));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% OTHER SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
disp(['Size of PSF ~= ' num2str(IMGSIZE_REF) ' [microlens pitch]' ]);
IMG_HALFWIDTH = max( Nnum*(IMGSIZE_REF + 1), 2*Nnum);
disp(['Size of IMAGE = ' num2str(IMG_HALFWIDTH*2*OSR+1) 'X' num2str(IMG_HALFWIDTH*2*OSR+1) '' ]);
x1space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; 
x2space = (pixelPitch/OSR)*[-IMG_HALFWIDTH*OSR:1:IMG_HALFWIDTH*OSR]; 
x1length = length(x1space);
x2length = length(x2space);

x1MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2];
x2MLspace = (pixelPitch/OSR)* [-(Nnum*OSR-1)/2 : 1 : (Nnum*OSR-1)/2];
x1MLdist = length(x1MLspace);
x2MLdist = length(x2MLspace);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% FIND NON-ZERO POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%
validpts = find(objspace>eqtol);
numpts = length(validpts);
[p1indALL p2indALL p3indALL] = ind2sub( size(objspace), validpts);
p1ALL = x1objspace(p1indALL)';
p2ALL = x2objspace(p2indALL)';
p3ALL = x3objspace(p3indALL)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%% DEFINE ML ARRAY %%%%%%%%%%%%%%%%%%%%%%%%% 
MLARRAY = calcML(fml, k0, x1MLspace, x2MLspace, x1space, x2space); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%% Alocate Memory for storing PSFs %%%%%%%%%%%   
LFpsfWAVE_STACK = zeros(x1length, x2length, numpts);
psfWAVE_STACK = zeros(x1length, x2length, numpts);
disp(['Start Calculating PSF...']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% PROJECTION FROM SINGLE POINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
centerPT = ceil(length(x1space)/2);
halfWidth =  Nnum*(IMGSIZE_REF + 0 )*OSR;
centerArea = (  max((centerPT - halfWidth),1)   :   min((centerPT + halfWidth),length(x1space))     );

disp(['Computing PSFs (1/3)']);
for eachpt=1:numpts,        
    p1 = p1ALL(eachpt);
    p2 = p2ALL(eachpt);
    p3 = p3ALL(eachpt);
    
    IMGSIZE_REF_IL = ceil(IMGSIZE_REF*( abs(p3)/p3max));
    halfWidth_IL =  max(Nnum*(IMGSIZE_REF_IL + 0 )*OSR, 2*Nnum*OSR);
    centerArea_IL = (  max((centerPT - halfWidth_IL),1)   :   min((centerPT + halfWidth_IL),length(x1space))     );
    disp(['size of center area = ' num2str(length(centerArea_IL)) 'X' num2str(length(centerArea_IL)) ]);    
    
    %% excute PSF computing funcion
    [psfWAVE LFpsfWAVE] = calcPSF(p1, p2, p3, fobj, NA, x1space, x2space, pixelPitch/OSR, lambda, MLARRAY, d, M, n,  centerArea_IL);
    psfWAVE_STACK(:,:,eachpt)  = psfWAVE;
    LFpsfWAVE_STACK(:,:,eachpt)= LFpsfWAVE;    

end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%% Compute Light Field PSFs (light field) %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
x1objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];
x2objspace = (pixelPitch/M)*[-floor(Nnum/2):1:floor(Nnum/2)];
XREF = ceil(length(x1objspace)/2);
YREF = ceil(length(x1objspace)/2);
CP = ( (centerPT-1)/OSR+1 - halfWidth/OSR :1: (centerPT-1)/OSR+1 + halfWidth/OSR  );
H = zeros( length(CP), length(CP), length(x1objspace), length(x2objspace), length(x3objspace) );

disp(['Computing LF PSFs (2/3)']);
for i=1:length(x1objspace)*length(x2objspace)*length(x3objspace) ,
    [a, b, c] = ind2sub([length(x1objspace) length(x2objspace) length(x3objspace)], i);  
     psfREF = psfWAVE_STACK(:,:,c);  
     
     psfSHIFT= im_shift2(psfREF, OSR*(a-XREF), OSR*(b-YREF) );
     [f1,dx1,x1]=fresnel2D(psfSHIFT.*MLARRAY, pixelPitch/OSR, d,lambda);
     f1= im_shift2(f1, -OSR*(a-XREF), -OSR*(b-YREF) );
     
     xmin =  max( centerPT  - halfWidth, 1);
     xmax =  min( centerPT  + halfWidth, size(f1,1) );
     ymin =  max( centerPT  - halfWidth, 1);
     ymax =  min( centerPT  + halfWidth, size(f1,2) );
 
     f1_AP = zeros(size(f1));
     f1_AP( (xmin:xmax), (ymin:ymax) ) = f1( (xmin:xmax), (ymin:ymax) );
     [f1_AP_resize, x1shift, x2shift] = pixelBinning(abs(f1_AP.^2), OSR);           
     f1_CP = f1_AP_resize( CP - x1shift, CP-x2shift );

     H(:,:,a,b,c) = f1_CP;

end

for aa=1:size(H,3)
    for multiWDF=1:size(H,4)
        for kk=1:size(H,5)
            temp=H(:,:,aa,multiWDF,kk);
            H(:,:,aa,multiWDF,kk)= H(:,:,aa,multiWDF,kk)./sum(temp(:));
        end
    end
end

x1space = (pixelPitch/1)*[-IMG_HALFWIDTH*1:1:IMG_HALFWIDTH*1];
x2space = (pixelPitch/1)*[-IMG_HALFWIDTH*1:1:IMG_HALFWIDTH*1]; 
x1space = x1space(CP);
x2space = x2space(CP);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%% Clear variables that are no longer necessary %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear LFpsfWAVE_STACK;
clear LFpsfWAVE_VIEW;
clear psfWAVE_VIEW;
clear LFpsfWAVE;
clear PSF_AP;
clear PSF_AP_resize;
clear PSF_CP;
clear f1;
clear f1_AP;
clear f1_AP_resize;
clear f1_CP;
clear psfREF;
clear psfSHIFT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tol = 0.005;
for i=1:size(H,5),
   H4Dslice = H(:,:,:,:,i);
   H4Dslice(find(H4Dslice< (tol*max(H4Dslice(:))) )) = 0;
   H(:,:,:,:,i) = H4Dslice;   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%% Estimate PSF size again  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
centerCP = ceil(length(CP)/2);
CAindex = zeros(length(x3objspace),2);
for i=1:length(x3objspace),
    IMGSIZE_REF_IL = ceil(IMGSIZE_REF*( abs(x3objspace(i))/p3max));
    halfWidth_IL =  max(Nnum*(IMGSIZE_REF_IL + 0 ), 2*Nnum);
    CAindex(i,1) = max( centerCP - halfWidth_IL , 1);
    CAindex(i,2) = min( centerCP + halfWidth_IL , size(H,1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%% Convert to phase-space PSF  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
IMGsize=size(H,1)-mod((size(H,1)-Nnum),2*Nnum);
psf=zeros(IMGsize,IMGsize,Nnum,Nnum,size(H,5)); %% phase-space PSF
for z=1:size(H,5)
 
    sLF=zeros(IMGsize,IMGsize,Nnum,Nnum); %% scanning light field images
    index1=round(size(H,1)/2)-fix(size(sLF,1)/2);
    index2=round(size(H,1)/2)+fix(size(sLF,1)/2);
    for ii=1:size(H,3)
        for jj=1:size(H,4)
            sLF(:,:,ii,jj)=im_shift3(squeeze(H(index1:index2,index1:index2,ii,jj,z)),ii-((Nnum+1)/2), jj-(Nnum+1)/2);
        end
    end
    
    multiWDF=zeros(Nnum,Nnum,size(sLF,1)/size(H,3),size(sLF,2)/size(H,4),Nnum,Nnum); %% multiplexed phase-space
    for i=1:size(H,3)
        for j=1:size(H,4)
            for a=1:size(sLF,1)/size(H,3)
                for b=1:size(sLF,2)/size(H,4)
                    multiWDF(i,j,a,b,:,:)=squeeze(  sLF(  (a-1)*Nnum+i,(b-1)*Nnum+j,:,:  )  );
                end
            end
        end
    end
    WDF=zeros(  size(sLF,1),size(sLF,2),Nnum,Nnum  ); %% multiplexed phase-space
    for a=1:size(sLF,1)/size(H,3)
        for c=1:Nnum
            x=Nnum*a+1-c;
            for b=1:size(sLF,2)/size(H,4)
                for d=1:Nnum
                    y=Nnum*b+1-d;
                    WDF(x,y,:,:)=squeeze(multiWDF(:,:,a,b,c,d));
                end
            end
        end
    end
    psf(:,:,:,:,z)=WDF;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp(['Saving PSF matrix file...']);
save('psf',['psf_M',num2str(M),'_NA',num2str(NA),'zmin',num2str(zmin),'u_zmax',num2str(zmax),'u.mat'],'-v7.3');
disp(['PSF computation complete.']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
