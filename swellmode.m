function DA = swellmode(Partic,D,target)
% K is the size increment for ridging/rafting
% We send our area into the gridbox of size fsize ~ k * lambda
% This is the GAIN. Due to other guys rafting in. The loss due to you
% rafting out is included before

numbins = length(D);
DA = 0*D; 

for i = 1:numbins
    %Find target
    %shift = (D - fsize);
    %[~,target] = min(shift>0);
    %Effectively performs the role of the delta function. Put the area into
    %the box that is smaller than the maximal cracked size. 
    
    %Put participation area in the target
    DA(target(i)) = DA(target(i)) + Partic(i); %No normalization required, I think. 
    
end

end
    

