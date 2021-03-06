function DA = redist_fsd(A,epsdot,f,D,shiftra,shiftri)

%This function controls the binned redistribution, following the work of
%Thorndike and others. Redistribution involves both ridging and rafting, as
%well as a dependence on the homogeneity of the thickness distribution. 

%The redistribution function looks like:
% Psi = |epsdot|*(A*delta(d) + B*(-P(r)a(r) + f*w_ra(r) +(1-f)*w_ri(r)))

% The w functions are predefined rafting and ridging distributors
% f is the homogeneity of the thickness distribution. 
% epsdot is the strain rate tensor
% A and B are the lead opening and lead closing coefficients, respectively
% P(r) is the participation function

% This code does not create open water explicity. If you would like this
% functionality, use Redist_OW_fsd.m which does this with variable
% open-water fractions k1 and k2.

mult = 10; %Artificial multiplier incase we want to boost redist

numbins = length(A);

inv = eig(epsdot);

% Strain Rate Invariants - Properly normalized
% Normalization doesn't affect the closing/opening ratio, thankfully

eps1 = (1/sqrt(2))*(inv(1) + inv(2)); %First strain rate invariant: the divergence
eps2 = (1/sqrt(2))*(inv(1) - inv(2)); %Second strain rate invariant: ~ the shear

strainmag = norm(epsdot,'fro'); %Frobenius Norm = root of sum of squares of entries

theta = atan(eps2/eps1); %Ratio of shear to divergence


%Opening and closing coefficients

%Just to begin with
leadopen = .5*(1 + cos(2*theta));
leadclose = .5*(1 - cos(2*theta));

% establish the participation function for each gridbox
% This operates on the idea that after a certain percentage, ice doesn't
% participate in rafting and ridging as (1 - f) where f = cum(A)/A*

cumA = cumsum(A); %This is the cdf


% ***This is the cdf-based participation function

Astar = .8; %From G* = .15 in Thorndike, area less likely to cut off

Partic1 = max((1 - cumA/Astar),0);
if sum(Partic1) > 1
Partic1 = Partic1/sum(Partic1);
end

% *** This participation function says we interact with our own class


%This would be normalized in a continuous setting
%Can't participate more than 1 time in an encounter
%It is not here, so we will just renormalize ourselves
% %Only if we get too much participation. Maybe? Not accounted for in
% Thorndike
% Should never be able to lose more mass than you have

Partic = Partic1.*A; %How much is lost in each cat. 

DA = 0*Partic; 

DA(1) = leadopen;

%Loss due to Participation
DAr = -Partic;

%DAr(1) = 0;
%% Ridging Mode
k = sqrt(5); % No open water creation explicit in redistribution. 

DAridge = (1-f)*redistmode(shiftri,k,D,Partic);% + .25*redistmode(min(shiftri+1,nbins),k,D,Partic) + .25*redistmode(max(shiftri-1,2),k,D,Partic));

%This has been smoothed. 
%% Rafting Mode

k = sqrt(3/2);

DAraft = f*redistmode(shiftra,k,D,Partic);

DA = mult*strainmag*leadclose*(DAr + DAridge + DAraft);

%% If we have no open water, redistribution must be adjusted so it stays that way

DA(1) = -sum(DA(2:numbins));

[C,I] = min(A + DA);


if C < 0 
    DAo = DA(1);
    DA(1) = DAo - C; 
    frac = abs(DA(1)/DAo);
    DA(2:numbins) = frac*DA(2:numbins);
end

 end

