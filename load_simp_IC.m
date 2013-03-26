function [A0,T,D,H,g,f,shiftra,shiftri,Yg,lambda,Spectrum,target,forcingperiod] = load_simp_IC(nbins)
 
D = logspace(0,2,nbins); %Bin Values, Include 0 and 1
 A0(1) = .5;
 for i = 2:nbins
     A0(i) = .5*A0(i-1);
 end
 A0(1) = .1; % No Open Water in Pack
 %Power Law w/ exponent
% A0(1) = .5; %Open water fraction
% A0(round(nbins/2)) = .25; %put a bunch in one part
 A0 = A0/sum(A0);

 T =-3; %Init Temp
 
forcingperiod = 500;
 
 
Have = 1; %m, average thickness
%Thickness Categories
H = logspace(0,log(5),nbins);

%Thickness distribution (not updated, as of now)

% Per Torge Martin, use a lognormal pdf with mean ln(h) - 1/32
% and deviation .25
mu = log(Have) - 1/32;
sigma = .25;
g = lognpdf(H,mu,sigma);
g = g/sum(g);

% Homogeneity
f = .5;

%This is just an estimate till I come up with a good way of defining the
%homogeneity of the ice pack.

% Peak Wave Period (s)
Tp = 10;

% Coefficients for rafting/ridging
% k1 = sqrt(3/2);
% k2 = sqrt(6/5);

k1 = sqrt(3/2);
k2 = sqrt(5);

%Shift matrices for rafting/ridging
    for i = 1:nbins
        %Find target
        rashift = (D  - k1*D(i)); %Shift due to rafting
        rishift = (D - k2*D(i)); %Shift due to ridging
    
        RA = rashift > 0;
        RI = rishift > 0;
        
        if max(RA) == 0
            shiftra(i) = nbins;
        else
        shiftra(i) = find(RA == 1, 1 );
        end
        
        if max(RI) == 0
            shiftri(i) = nbins;
        else
        shiftri(i) = find(RI == 1, 1 );
        end
     % [~,shiftra(i)] = min(rashift>0); %Effectively performs the role of the delta function
     % [~,shiftri(i)] = min(rishift>0);
    end
    
 %UNCOMMENT THIS LINE IF WE WANT THE OPEN WATER NOT TO RIDGE/RAFT
    %shiftri(1) = 1;
%shiftra(1) = 1;
%% Doing stuff for the swell fracture

nbins = length(D); %# of area bins
grav = 9.81; %m/s^2
Hs = 0; 
epsc = 5e-5;
sigc = .8e6;
rhobar = .5*(1025 + 922.5);

%% Incident Wave Spectrum

%Twice as many possible wavelengths as bins
%This analysis only really works for periods between 6 and 16 seconds as
%outlined by Kohout and Meylan. But we will extend it somwhat down to a
%wavelength of 10 m.

lambda = linspace(2*D(1) + .00001,2*D(nbins) + .00001,2*nbins);

%Deep water, ice-free dispersion relation
per = sqrt(2*pi*lambda/grav); 

%Require Significant wave height in BS
constBS = 1.25*Hs^2/(8*pi*Tp^4);
%No significant wave height in PMS
constPMS = (8.1e-3)*(grav^2)*(2*pi)^(-5);

const = constPMS; % Use PMS for now

Spectrum = const*per.^(5).*(exp(-1.5*(per/Tp).^4)); %May need to normalize

Amp = sqrt(4*pi*Spectrum./per); %Amplitude at each period

Sp = Spectrum/sum(Spectrum); %Normalized spectrum pdf

%% Calculate Critical Amplitude for the thickness distribution

for i = 1:length(H)
    %For each thickness, 
    Aecrit(i,:) = epsc*lambda.^2./(2*pi*H(i)); 
    Ascrit(i,:) = 2*pi*H(i)^2*sigc./(3*grav*rhobar*lambda.^2);
    % Find minimal amplitude to crack the ice
end

Acrit = min(Aecrit,Ascrit);
    
    
%% Calculate Stress/Strain Failure Percentage from H
Yg = 0*lambda;
% At each lambda and h, does it meet the critical threshold?

sig = bsxfun(@gt,Amp,Acrit);

% If so, add to Yg the value of its percentage of the thickness dist

Yg = g*sig; % THIS MAY BE WRONG!!! CHECK 

%Now we have Yg(lambda)

for i = 1:length(lambda)
        %Find target
        shift = (D - .5*lambda(i));
        SH = shift < 0;
        target(i) = find(SH == 1,1,'last');
end

end