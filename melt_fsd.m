 function Melt = melt_fsd(A,T,D)
% function Melt = melt_fsd(A,T)
% This function computes the rate of lateral melting of ice floes based on
% the temperature difference between the ice and the ocean. 

%% Specify Constants

Lf = 334000; %J/kg
rhoi = 900; %kg/m^3
ki = 2.03e5; %Thermal Conductivity J/(degC S meter)
Tfr = -1.96; %Temp of Freezing (deg C)

eoff = .05; %Renormalization Offset

%% Compute Thermodynamic Melt Rate

nbins = length(D); 
open = A(1); %Open Water

Cond = open.*(ki*(T - Tfr))/(Lf*rhoi);

% Conductive growth rate,
% Limiting when there is no open water

% Thermodynamic Areal Growth Rate

% This so that each slice gets an equal amount of areal growth

Therm = Cond + 0*D; 

%.*D./(pi*(D+eoff).^2);

%% Compute Fluxes 
dd(1) = 0;
dfa(1) = 0;
for i = 2:nbins-1 
    %No thermo growth of open water... open water is just fixed to be the opposite of the other growth
    
    % delta d
    dd(i) = D(i+1) - D(i);
   
    % "Advective" flux in diameter space
    dfa(i) = A(i+1).*Therm(i+1) - A(i).*Therm(i);
end

% Handle the final bin
dd(nbins) = dd(nbins - 1);
dfa(nbins) = A(nbins).*Therm(nbins); 
    
    
Melt = dfa./dd;
Melt(1) = -sum(Melt(2:nbins));
%Done!    
    


    
