function [] = imwriteTFSK(Volume,name)
%% write 3D image file 
%% supported formats: single/double/uint8/uint16
%% Input:
% @Volume: the 3D image  
% @name: the file path
%
%   [1]  JIAMIN WU, ZHI LU and DONG JIANG etc,
%        Iterative tomography with digital adaptive optics permits hour-long intravital observation of 3D subcellular dynamics at millisecond scale
%        Cell, 2021. 
%
%    Contact: ZHI LU (luz18@mails.tsinghua.edu.cn)
%    Date  : 10/24/2020

t = Tiff(name,'w'); % Filename by variable name
tagstruct.ImageLength       = size(Volume,1);
tagstruct.ImageWidth        = size(Volume,2);
tagstruct.Photometric       = Tiff.Photometric.MinIsBlack;
if(class(Volume) == 'single')
    tagstruct.BitsPerSample	= 32;
    tagstruct.SampleFormat	= Tiff.SampleFormat.IEEEFP;
elseif(class(Volume) == 'double')
    warning('ImageJ may not support double/64-bit tiff!');
    tagstruct.BitsPerSample	= 64;
    tagstruct.SampleFormat	= Tiff.SampleFormat.IEEEFP;
elseif(class(Volume) == 'uint8')
    tagstruct.BitsPerSample	= 8;
    tagstruct.SampleFormat	= Tiff.SampleFormat.UInt;
elseif(class(Volume) == 'uint16')
    tagstruct.BitsPerSample	= 16;
     tagstruct.SampleFormat	= Tiff.SampleFormat.UInt;
end
tagstruct.SamplesPerPixel	= 1;
tagstruct.RowsPerStrip      = 16;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software          = 'MATLAB';
setTag(t,tagstruct)
write(t,squeeze(Volume(:,:,1)));
for i=2:size(Volume,3) % Write image data to the file
    writeDirectory(t);
    setTag(t,tagstruct)
    write(t,squeeze(Volume(:,:,i))); % Append
end
close(t);
end

