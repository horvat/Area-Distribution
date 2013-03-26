function [A0,T,D,H,g,f,shiftra,shiftri,epsri,epsra,Yg,lambda,Spectrum,target,forcingperiod,Multri,Multra] = load_simp_IC(nbins,version)
 
D = linspace(1,200,nbins); %Bin Values, Include 0 and 1

Dquart = D(ceil(nbins/3)); 
Dquart2 = D(ceil(nbins/2)); 

A0 = normpdf(D,Dquart,Dquart/10);% + normpdf(D,Dquart2,Dquart/10); 
 
A0 = A0/sum(A0); 
 
 
forcingperiod = 500;

T =-3; %Init Temp
 
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
f = 0;

%This is just an estimate till I come up with a good way of defining the
%homogeneity of the ice pack.

% Peak Wave Period (s)
Tp = 10;

% Coefficients for rafting/ridging
% k1 = sqrt(3/2);
% k2 = sqrt(6/5);

if version == 2
    
    % These are now matrices since this is a two-floe interaction
    
    shiftra = zeros(nbins,nbins);
    shiftri = shiftra;
    
    % The assumption is: two pieces come together
    % The smallest of the two has some part of it turned into the ridge
    % (This is no different than saying both ridge, since they'll combine)
    % r percent of the original piece ridges/rafts
    % This turns into a ridge of height k times the original height
    % Therefore the original piece has an area which is now: 
    % A1 = (1-r)A1 + rA1 --> (1-r)A1 + (r/k)A1 = (1 - (k-1)r/k)A1 = A1'
    % And the total area becomes A2+A1 --> A2 + A1'
    % Plus an open water term
    
 % In ridging, k ~ 5 and r ~ .25 So a quarter of the incident piece is
 % turned into a ridge
 % In rafting, k ~ 2 and r ~ 1 So the entire piece is turned into a double
 % thick piece. Or, the whole piece uplifts onto the other piece and half
 % of it remains. Same idea. 
 
    % This gives for the new diameter:
    % dnew = (d1^2 + (1 - (k-1)r/k)d2^2)^(1/2)
    kri = 5;
    kra = 2; 
    rri = .25;
    rra = 1;
    
    fracri = 4/5;
    fracra = 1/2;
    
    epsri = sqrt(1 - fracri);
    epsra = sqrt(1 - fracra); 
    for i = 2:nbins
        for j = i:nbins
            dnewri = sqrt(D(j)^2 + fracri*D(i)^2); 
            dnewra = sqrt(D(j)^2 + fracra*D(i)^2);
            
            
            rashift = (D - dnewra); 
            rishift = (D - dnewri); 
            
            RA = rashift > 0; 
            RI = rishift > 0;
            
            if max(RA) == 0; 
                shiftra(i,j) = nbins;
            else
                shiftra(i,j) = find(RA == 1, 1 );
            end
            
            if max(RI) == 0
                shiftri(i,j) = nbins;
            else
                shiftri(i,j) = find(RI == 1, 1 );
            end
        end
    end
 
else
    
    % Original Version shift matrices
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
    
end
    
 %UNCOMMENT THIS LINE IF WE WANT THE OPEN WATER NOT TO RIDGE/RAFT
    %shiftri(1) = 1;
%shiftra(1) = 1;

Multri = 1./(bsxfun(@plus,D.^2,(D').^2)./bsxfun(@plus,fracri*D.^2,(D').^2)); 
Multra = 1./(bsxfun(@plus,D.^2,(D').^2)./bsxfun(@plus,fracra*D.^2,(D').^2)); 


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