function znew = Timeweighted_interp(ref,index,a1,a2,Nshift)
%% Input:
% @ref: the raw phase-space image  
% @index: the time-weighted distribution
% @a1: the scanning index in the first dimension
% @a2: the scanning index in the second dimension
% @Nshift: the sampling points of a single scanning period
%% Output:
% new_im: the shifted 3D image  
%
%    [1]  JIAMIN WU, ZHI LU, DONG JIANG and YUDUO GUO.etc,
%         3D observation of large-scale subcellular dynamics in vivo at the millisecond scale
%         in BioRxiv, 2019. 
%
%    Contact: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

% distanced weighting
index_ref =index;
index = index_ref./sum(sum(index_ref(1:3,1:3)));
stack = zeros(size(ref,1),size(ref,2),Nshift^2);
for i = 1:Nshift^2
     tmp = index(a1(i),a2(i))*imresize(ref(a1(i):Nshift:end,a2(i):Nshift:end),[size(ref,1),size(ref,2)],'cubic');
     stack(:,:,i) = im_shift3(tmp,(-round(Nshift/2)+a1(i)),(-round(Nshift/2)+a2(i)));
end
z_inter = sum(stack,3);

% generate timeweighted phase-space
znew = (1-index_ref).*z_inter + index_ref .* ref;
end














