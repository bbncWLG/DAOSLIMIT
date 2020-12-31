function Xguess = deconvRL(maxIter, Xguess,WDF, psf,weight,DAO,GPUcompute)
% Deconvolution for the single-frame multiplexed phase-space data
%% Input:
% @maxIter: the maximum iteration number 
% @Xguess: the estimated volume
% @WDF: multiplexed phase-space data
% @psf: phase-space PSF
% @weight: the weight numbers corresponding to the all spatial frequencies
% @AO: digital adaptive optics on/off
% @GPUcompute: GPU on/off
%% Ouput: 
% output: the high resolution 3D volume by deconvolution  
%
%   [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%        3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%        in BioRxiv, 2019. 

% Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
% Date  : 10/24/2020

Nnum=size(psf,3); %% number of virtual pixels
if GPUcompute==1  %% GPU on
    WDF = gpuArray(single(WDF));
    weight=gpuArray(single(weight));
    Xguess=gpuArray(single(Xguess));
end

% 3D deconvolution order
index1=[1,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,4,4,4,4,4,4,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,5,5,5,5,5,6,7,8,9,9,9,9,9,8,7,6,6,6,6,7,8,8,8,7,7];
index2=[1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,13,13,13,13,13,13,13,13,13,13,12,11,10,9,8,7,6,5,4,3,2,1,1,1,1,1,1,1,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,12,12,12,12,12,12,12,12,12,12,11,10,9,8,7,6,5,4,3,2,2,2,2,2,2,2,2,2,2,3,4,5,6,7,8,9,10,11,11,11,11,11,11,11,11,11,10,9,8,7,6,5,4,3,3,3,3,3,3,3,3,4,5,6,7,8,9,10,10,10,10,10,10,10,9,8,7,6,5,4,4,4,4,4,4,5,6,7,8,9,9,9,9,9,8,7,6,5,5,5,5,6,7,8,8,8,7,6,6,7];

% iteractive DAO deconvolution algorithm
for iter=1:maxIter
    tic;
    % Digital adaptive optics to estimate aberration
    if DAO>0 % DAO on
        sidelobe=round(0.04*size(WDF,1)); %% reserved image border
        Nb=1;
        map_wavshape=zeros(Nnum,Nnum,Nb,Nb,2);     
        Sb=fix( 0.9*size(WDF,1)/(Nb)/2 )*2+1;
        border=(size(WDF,1)-Sb*Nb)/2;
        weight_mask=single(im2bw(gather(weight),1e-5));
        [Sx,Sy]=meshgrid(-fix(Nnum/2):fix(Nnum/2),-fix(Nnum/2):fix(Nnum/2));
        if iter>1
            for i=1:169
                u=index1(i);
                v=index2(i);
                if weight(u,v)==0
                    continue;
                else
                    if GPUcompute==1
                        psf_uv=gpuArray(single(squeeze(psf(:,:,u,v,:))));
                        forwardFUN = @(Xguess) forwardProjectGPU( psf_uv, Xguess );
                    else
                        psf_uv=single(squeeze(psf(:,:,u,v,:)));
                        forwardFUN = @(Xguess) forwardProjectACC( psf_uv, Xguess );
                    end
                    HXguess=forwardFUN(Xguess); 
                    
                    % spatial partition to estimate wavefronts
                    for uu=1:Nb
                        for vv=1:Nb
                            % tiled correlation with partitioning
                            sub_HXguess=HXguess(border+(uu-1)*Sb+1:border+uu*Sb,border+(vv-1)*Sb+1:border+vv*Sb);
                            sub_WDF=WDF(border+(uu-1)*Sb+1-sidelobe:border+uu*Sb+sidelobe,border+(vv-1)*Sb+1-sidelobe:border+vv*Sb+sidelobe,u,v);
                            corr_map=gather(normxcorr2(sub_HXguess,sub_WDF));
                            [testa,testb]=find(corr_map==max(corr_map(:)));
                            map_wavshape(u,v,uu,vv,1)=testa-size(sub_WDF,1)+sidelobe;
                            map_wavshape(u,v,uu,vv,2)=testb-size(sub_WDF,2)+sidelobe;       
                        end
                    end
                    for uu=1:Nb
                        for vv=1:Nb
                            cx=map_wavshape(7,7,uu,vv,1);
                            cy=map_wavshape(7,7,uu,vv,2);
                            map_wavshape(:,:,uu,vv,1)=(squeeze(map_wavshape(:,:,uu,vv,1))-cx).*weight_mask;
                            map_wavshape(:,:,uu,vv,2)=(squeeze(map_wavshape(:,:,uu,vv,2))-cy).*weight_mask;
                        end
                    end
                    for uu=1:Nb
                        for vv=1:Nb
                            for u=1:Nnum
                                for v=1:Nnum
                                    map_wavshape(u,v,uu,vv,1)=min(max(map_wavshape(u,v,uu,vv,1),-sidelobe),sidelobe);%% limit the shifted range
                                    map_wavshape(u,v,uu,vv,2)=min(max(map_wavshape(u,v,uu,vv,2),-sidelobe),sidelobe);
                                end
                            end
                        end
                    end
                    % remove Zernike defocus items from wavefronts
                    for uu=1:Nb
                        for vv=1:Nb
                            k1 = Sy.*squeeze(map_wavshape(:,:,uu,vv,1))+Sx.*squeeze(map_wavshape(:,:,uu,vv,2));
                            k2 = Sx.*Sx+Sy.*Sy;
                            k=sum(k1(:))/sum(k2(:));
                            map_wavshape(:,:,uu,vv,1)=squeeze(map_wavshape(:,:,uu,vv,1))-k*Sy;
                            map_wavshape(:,:,uu,vv,2)=squeeze(map_wavshape(:,:,uu,vv,2))-k*Sx;
                            for u=1:Nnum
                                for v=1:Nnum
                                    map_wavshape(u,v,uu,vv,1)=min(max(map_wavshape(u,v,uu,vv,1),-sidelobe),sidelobe);
                                    map_wavshape(u,v,uu,vv,2)=min(max(map_wavshape(u,v,uu,vv,2),-sidelobe),sidelobe);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    % Volume update
    for i=1:169        
        u=index1(i);
        v=index2(i);       
        if weight(u,v)==0
            continue;
        else   
            if DAO>0 && iter>1
                % correct aberration
                map_wavshape_x=squeeze(map_wavshape(u,v,:,:,1));
                map_wavshape_y=squeeze(map_wavshape(u,v,:,:,2));
                map_wavshape_x1=imresize(map_wavshape_x,[round(size(WDF,1)/3),round(size(WDF,2)/3)],'nearest');
                map_wavshape_y1=imresize(map_wavshape_y,[round(size(WDF,1)/3),round(size(WDF,2)/3)],'nearest');
                map_wavshape_xx=imresize(map_wavshape_x1,[size(WDF,1),size(WDF,2)],'cubic');
                map_wavshape_yy=imresize(map_wavshape_y1,[size(WDF,1),size(WDF,2)],'cubic');
                [coordinate1,coordinate2]=meshgrid(1:size(WDF,1),1:size(WDF,2));
                if GPUcompute==1
                    WDF_uv=gpuArray(interp2(coordinate1,coordinate2,WDF(:,:,u,v),coordinate1+map_wavshape_yy,coordinate2+map_wavshape_xx,'cubic',0));
                else
                    WDF_uv=interp2(coordinate1,coordinate2,WDF(:,:,u,v),coordinate1+map_wavshape_yy,coordinate2+map_wavshape_xx,'cubic',0);
                end
            else
                WDF_uv=WDF(:,:,u,v);
            end
            if GPUcompute==1
                 psf_uv=gpuArray(single(squeeze(psf(:,:,u,v,:))));
                 forwardFUN = @(Xguess) forwardProjectGPU( psf_uv, Xguess );
                 backwardFUN = @(projection) backwardProjectGPU(psf_uv, projection );
                 uniform_matrix = gpuArray(single (ones(  size(WDF,1),size(WDF,2) )));
            else
                psf_uv=single(squeeze(psf(:,:,u,v,:)));
                forwardFUN = @(Xguess) forwardProjectACC( psf_uv, Xguess );
                backwardFUN = @(projection) backwardProjectACC(psf_uv, projection );
                uniform_matrix = single (ones(  size(WDF,1),size(WDF,2) ));
            end 
            % RL deconvolution
            HXguess=forwardFUN(Xguess);      
            errorEM=squeeze(WDF_uv)./HXguess;
            errorEM(~isfinite(errorEM))=0;
            XguessCor = backwardFUN(errorEM) ; 
            Htf=backwardFUN( uniform_matrix );
            Xguess_add=Xguess.*XguessCor./Htf;
            clear Htf;clear XguessCor;
            Xguess_add(find(isnan(Xguess_add))) = 0;
            Xguess_add(find(isinf(Xguess_add))) = 0;
            Xguess_add(Xguess_add<0 ) = 0;
            Xguess=Xguess_add.*weight(u,v)+(1-weight(u,v)).*Xguess;
            clear Xguess_add;
            Xguess(find(isnan(Xguess))) = 0;
            Xguess(Xguess<1e-4)=0;
        end
    end
    ttime = toc;
    disp(['  iter ' num2str(iter) ' | ' num2str(maxIter),', took ' num2str(ttime) ' secs']);    
%     imwriteTFSK(gather(Xguess),['AO1_2/iter',num2str(iter),'.tif']);
end