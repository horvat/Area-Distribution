%function [D A DA normDA] = areadistribution 
% Version 0.3


%%

% This script is the driver for the area distribution code. It simulates
% the Floe Size Distribution as it is outlined in Chris Horvat's Proposal. 
% It includes the effects of swell fracture, lateral melting, rafting and
% ridging as well as advection. The model is binned and dynamics are on the
% Arakawa A-grid. 

% The thickness distribution is assumed to be given as log-normal (see
% Torge Martin's work with MITGCM) with a log standard deviation of .24 and
% a log mean of ln(h) - 1/32, where h is the model-given mean thickness

% Sea ice dynamics are assumed to be donated by some benevolent GCM. The
% stress tensor employs the viscous-plastic constitutive law of
% Hibler(1980). 

%The area distribution feeds back into the coupled system via its influence
%in air-sea fluxes. This is not included yet, the purpose of this model is
%to diagnose the floe size distribution. 

%The Model equations are as foll 
% Advection + Melt = Redistribution + Swell Fracture

%At first, the model will be examined in one grid cell. Boundary conditions
%are: Pack ice at the top of the cell, open water at the bottom. 
clear all
%close all

%% Set parameters
Ttotal = 500000; %Model time
dt = 1; %Time interval
dumpfreq = dt; %Frequency of saving
nt = Ttotal/dt; %Number of iterations
Time = linspace(0,dt*nt,nt); %Time vector
nbins = 60; %Number of bins for Area Distribution
tol = 1e-6;

disp(sprintf('Computing the FSD with a steady-state tolerance of %d for up to %d seconds',tol,Ttotal))

%% Initialization

% A = zeros(1,nbins); %Area Distribution. Just one gridbox for now. 

%[A(1,:),H,Vel,Hs,Tpeak,T,Vair,Voce,Epsdot,Sigma,D] = load_initial_conditions(nbins);
version = 2; 
% Describes the rafting/ridging version

[A0,T,D,H,g,f,shiftra,shiftri,epsri,epsra,Yg,lambda,Spectrum,target,forcingperiod,Multri,Multra] = load_simp_IC(nbins,version);
% We load in the rafting and ridging shifts because doing this in each
% redistribution step kills us timewise in execution

L = D(end); 

A = A0;
% A(1,:) is the initial area distribution
% H is the initial mean thickness of the ice field
% V is the initial ice velocity field
% H_s is the "significant wave height" for the Pierson-Moskovitz Spectrum
% Tpeak is the peak period for the P-M Spectrum
% T(1:nt) is the atmos. temperature field
% Vair(1:nt) is the atmos. velocity field
% Voce(1:nt) is the ocean velocity field
% Epsdot is the initial strain rate tensor
% Sigma is the initial Stress Tensor
% D is a vector giving the middle value for the diameter in each bin,
% except for D(1) = 0. 
% shiftri and shiftra describe, for a given initial piece, what it will
% raft or ridge into. they are matrices in the second version, just vectors
% in the first version. This is explained in load_simp_IC and redist_fsd_v2

[A1,A2,V1,V2,epsdot] = load_bc(D,0);
%A1,V1 is the area dist,velocity at the top (pack bc)
%A2 is the area dist,velocity at the bottom (ocean bc)

%% Grid
% The grid looks like
%
%      _________
%     |  PACK   |
%     |---------|
%  y  |   MIZ   |
%     |---------|
%     |  OCEAN  |
%     |_________|
%          
%          x



%% Updating Function
DA = 0;
stats = zeros(round(nt/dumpfreq),7);
sumstats = zeros(round(nt/dumpfreq),6);
sumstats2 = zeros(round(nt/dumpfreq),nbins);
Vals = zeros(round(nt/dumpfreq),nbins*4);
Vals = reshape(Vals,[round(nt/dumpfreq) 4 nbins]);
Save = zeros(4,nbins);
Save2 = zeros(1,nbins);
Aint = 0;
Melt = 0*D;
Redist = 0*D;
Advect = 0*D;
Swell = 0*D;
DA = ones(size(D));
t = 2;
tmax = 25000;
%A0 = A1;

while t < nt+1 && (norm(DA) > tol)
    [A1,A2,V1,V2,epsdot] = load_bc(D,0);
    DA = 0*D;
    % Compute and updat the ice velocity and stress/strain tensors
    %% [Sigma,Epsdot,Vel,H] = do_ice_dynamics(t,H,Vel,Sigma,Epsdot,Voce,Vair)   

    %% Advect Ice - Coded completely, no issues
    %Advect = .5*x_advect_fsd(A,D,V1,V2,A1,A2);       %.5 just because    
    Advect = 0*D;
    Aloss = sum(Advect);
    Advect(1) = Advect(1) - Aloss; 
    
    DA = DA + Advect; 
    
    
    
    %% Thermodynamic Melting - Coded
    %Melt = melt_fsd(A,T,D);
    Melt = 0*D;
    
    DA = DA + Melt; 
    
    
    
    
    % Mechanical Redistribution - Coded completely, no issues
    if version == 2 % More involved redistribution
            Redist = redist_fsd_v2(A,epsdot,f,D,shiftra,shiftri,epsri,epsra,Multri,Multra); 
           % Redist = horvsmooth(Redist1,A); 
    else
             Redist = redist_fsd(A,epsdot,f,D,shiftra,shiftri);
    end
    
    Redist(1) = (Redist(1) - Aloss); % This handles the divergence
   % Redist = 0*D; 
        
    DA = DA + Redist; 
    
    %% Swell Fracture - Coded, Runs, maybe not ready
    %Swell = swellfrac_fsd(A,D,lambda,Yg,target,Spectrum,dt);
    Swell = 0*D; 
    
    DA = DA + Swell; 

    %% We're Done!


     Anew = A + dt*DA;

    %% Compute relative magnitude of terms, stats, and bug control
    
    if min(Anew) < 0 
        disp(['Negative value at t = ' sprintf('%d',dt*t)])
        Ttotal = t;

        break
    else
        if max(isnan(Anew)~= 0)
            disp(['NaN value at t = ' sprintf('%d',dt*t)])
            break
        end
    end
    
    A = Anew; 
    
    Save = Save + dt*[Advect; Melt; Redist; Swell];
    
    Save2 = Save2 + dt*Swell; 
    
    if mod(t,dumpfreq) == 0
        stats(t/dumpfreq,:) = [norm(Advect) norm(Melt) norm(Melt)/norm(Advect) ...
         norm(Redist) norm(Redist)/norm(Advect) norm(Swell) ...
         norm(Swell)/norm(Advect)];
     
        sumstats(t/dumpfreq,:) = [A(1) sum(Advect) sum(Melt) ...
         sum(Redist) sum(Swell) sum(DA(2:end))];
     sumstats2(t/dumpfreq,:) = Swell';
     
    Vals(t/dumpfreq,:,:) = [Melt; Redist; Advect;Swell];

    end  
    t = t + 1;
end

%% Are we at a steady state or have we run out of time?

if norm(DA) < tol
    disp(sprintf('Reached Steady State by t = %d',t*dt));
else
    disp('Reached time limit')
end

t = round((t-1)/dumpfreq);

stats = stats(1:t,:);
Vals = Vals(1:t,:,:);

%% Plot everything

plotall = 1;

if plotall == 1
    
Amin = 1/(pi*D(nbins)^2); %Maximal Area

saveplots = 0;

Saver = Save./(repmat(A0,[4 1]));

Timer = linspace(0,dt*t,length(stats(1:t,1)));

plotFSD(D,[A+.0001;A0+.0001;A1+.0001],Timer', ... 
    [0*Timer + 1;stats(:,3)';stats(:,5)';stats(:,7)'], ... 
    [stats(:,1)';stats(:,2)';stats(:,4)';stats(:,6)'], ... 
    abs(Saver)/(dt*nt) + .0000001,squeeze(Vals(end,:,:)),t,saveplots);

%Tquit = round((t-1)/dumpfreq);

%plotbalance(D,[A;A0;A1],squeeze(Vals(end,:,:)),t);
 
end

%%
figure

plot(D,A-A0,'-ok','Linewidth',2)
hold on
plot(D,(norm(A-A0)/norm(A))*A,'-r','Linewidth',3); 
plot(D,(norm(A-A0)/norm(A))*A0,'-b','Linewidth',1); 


