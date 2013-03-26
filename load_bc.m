function [A1,A2,V1,V2,epsdot] = load_bc(D,t)
nbins = length(D);
A2 = ones(1,nbins)/(nbins); 
A1 = D.^(-.5) + D.^(-2)*sin(t);
A1(1) = 1;
A1 = A1/sum(A1);
A2 = A1;
V1 = -.1;% +.075*sin(t); %m/s down
V2 = -.1;%+ .075*cos(t); %m/s down
eps = (V1.^2 + V2.^2)^(1/2)/max(D);
exx = (V1 - V2)/max(D); %Vertical Divergence

epsdot = [exx .1*(exx+eps); .1*(exx+eps) .5*exx];
% Initial Strain Rate tensor. eps is used so that there will be shear even
% with no divergence. It will be small, though. 

end
