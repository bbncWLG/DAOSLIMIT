function MLARRAY = calcML(fml, k, x1MLspace, x2MLspace, x1space, x2space)
%% Robert's code
%% Reference:  Robert Prevede, Young-Gyu Yoon, Maximilian Hoffmann, Nikita Pak.etc. 
%% "Simultaneous whole-animal 3D imaging of neuronal activity using light-field microscopy " 
%% in Nature Methods VOL.11 NO.7|July 2014.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1length = length(x1space);
x2length = length(x2space);
x1MLdist = length(x1MLspace);
x2MLdist = length(x2MLspace);
x1center = find(x1space==0);
x2center = find(x2space==0);
x1centerALL = [  (x1center: -x1MLdist:1)  (x1center + x1MLdist: x1MLdist :x1length)];
x1centerALL = sort(x1centerALL);
x2centerALL = [  (x2center: -x2MLdist:1)  (x2center + x2MLdist: x2MLdist :x2length)];
x2centerALL = sort(x2centerALL);

zeroline = zeros(1, length(x2space) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
patternML = zeros( length(x1MLspace), length(x2MLspace) );
patternMLcp = zeros( length(x1MLspace), length(x2MLspace) );
for a=1:length(x1MLspace),
    for b=1:length(x2MLspace),        
        x1 = x1MLspace(a);
        x2 = x2MLspace(b);
        xL2norm = x1^2 + x2^2;
        

        patternML(a,b) = exp(-i*k/(2*fml)*xL2norm);   
        patternMLcp(a,b) = exp(-0.05*i*k/(2*fml)*xL2norm);  
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MLspace = zeros( length(x1space), length(x2space) );
MLcenters = MLspace;
for a=1:length(x1centerALL),
    for b=1:length(x2centerALL),
        MLcenters( x1centerALL(a), x2centerALL(b)) = 1;
    end
end
MLARRAY = conv2(MLcenters, patternML, 'same');
MLARRAYcp = conv2(MLcenters, patternMLcp, 'same');

MLARRAYcpANG = angle(MLARRAYcp);
MLARRAYcpANG = MLARRAYcpANG - min(min(MLARRAYcpANG)) + 0.0;
MLARRAYcpANGnorm = MLARRAYcpANG/max(max(MLARRAYcpANG));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%