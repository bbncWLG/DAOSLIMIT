function new_im = im_shift3(img, SHIFTX, SHIFTY)
%% Input:
% @img: the 3D image  
% @SHIFTX: the shift in the first dimension
% @SHIFTY: the shift in the second dimension
%% Output:
% new_im: the shifted 3D image  
%
%    Author: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020
%
%% modified from Robert's code (im_shift2.m)
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.

eqtol = 1e-10;

xlength = size(img,1);
ylength = size(img,2);

if abs(mod(SHIFTX,1))>eqtol | abs(mod(SHIFTY,1))>eqtol
   error('SHIFTX and SHIFTY should be integer numbers');
end

% if SHIFTX >= xlength | SHIFTY >= ylength,
%    error('SHIFTX  and SHIFTY should be smaller than size(img,1) and size(img,2), respectively');
% end


    
SHIFTX = round(SHIFTX);
SHIFTY = round(SHIFTY);

new_im = zeros(xlength, ylength,size(img,3));

if SHIFTX >=0 & SHIFTY >= 0,
    new_im( (1+SHIFTX:end), (1+SHIFTY:end) ,:) = img( (1:end-SHIFTX), (1:end-SHIFTY) ,:);
elseif SHIFTX >=0 & SHIFTY < 0,
    new_im( (1+SHIFTX:end), (1:end+SHIFTY) ,:) = img( (1:end-SHIFTX), (-SHIFTY+1:end) ,:);
elseif SHIFTX <0 & SHIFTY >= 0,
    new_im( (1:end+SHIFTX), (1+SHIFTY:end) ,:) = img( (-SHIFTX+1:end), (1:end-SHIFTY) ,:);
else
    new_im( (1:end+SHIFTX), (1:end+SHIFTY) ,:) = img( (-SHIFTX+1:end), (-SHIFTY+1:end) ,:);
end


end