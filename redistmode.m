function DA = redistmode(shift,k,D,Partic)
% K is the size increment for ridging/rafting
% We send our area into the spots specified by shift, which is like sending
% stuff to a multiple of k. 
% This is the GAIN. Due to other guys rafting in. The loss due to you
% rafting out is included before

%%
numbins = length(D);
DA = 0*D; 

% Open Water cannot distribute, so we start from 2
for i = 1:numbins
    
    %Find target
    target = shift(i);
    %Put participation area in the target
    DA(target) = DA(target) + Partic(i);
    
end
if max(DA(2:numbins)) ~= 0
s1 = sum(DA(2:numbins));
DA(2:numbins) = smooth(DA(2:numbins))';
s2 = sum(DA(2:numbins));
DA(2:numbins) = (s1/s2)*DA(2:numbins);
end

%end
    

