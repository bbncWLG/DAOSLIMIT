function Backprojection = backwardProject(H, projection, Nnum )
%% Robert's code
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.

x3length = size(H,5);
Backprojection = zeros(size(projection, 1), size(projection, 1), x3length );
for aa=1:Nnum,
    for bb=1:Nnum,
        for cc=1:x3length,
                     
            Ht = imrotate( H(:,:,aa,bb,cc), 180);            
            tempSlice = conv2(projection, Ht, 'same');            
            Backprojection((aa:Nnum:end) , (bb:Nnum:end),cc) = Backprojection((aa:Nnum:end) , (bb:Nnum:end),cc) + tempSlice( (aa:Nnum:end) , (bb:Nnum:end) );          
        end
    end
end

