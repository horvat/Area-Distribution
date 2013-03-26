%function [D A DA normDA] = areadistribution 
% Version 0.2

%% Changes from version 0.1

% Plotting
% - Added plotbalance.m, which plots the long-term balance of terms

% Redistribution
% - Changed k in redist_fsd to allow for open water formation
% - Added an artificial control in redistribution to handle A(1) = 0; 
% - Put a multiplier for redistribution inside the redist_f stuff. 

% Swell
% - Added smoothing of output

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
close all

%% Set parameters
Ttotal = 5000; %Model time
dumpfreq = 1; %Frequency of saving
dt = 1; %Time interval
nt = Ttotal/dt; %Number of iterations
Time = linspace(0,dt*nt,nt); %Time vector
nbins = 20; %Number of bins for Area Distribution
tol = 1e-5;

disp(sprintf('Computing the FSD with a steady-state tolerance of %d for up to %d seconds',tol,Ttotal))

%% Initialization

% A = zeros(1,nbins); %Area Distribution. Just one gridbox for now. 

%[A(1,:),H,Vel,Hs,Tpeak,T,Vair,Voce,Epsdot,Sigma,D] = load_initial_conditions(nbins);

[A0,T,D,H,g,f,shiftra,shiftri,Yg,lambda,Spectrum,target,forcingperiod] = load_simp_IC(nbins);
% We load in the rafting and ridging shifts because doing this in each
% redistribution step kills us timewise in execution


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
Vals = reshape(Vals,[round(nt/dumpfreq) 4 20]);
Save = zeros(4,nbins);
Save2 = zeros(1,nbins);


Temps = [linspace(-5,5,10000) linspace(5,-5,10000)];
for i = 1:length(Temps)
    T = Temps(i); 
    Aint = 0;
    Melt = 0*D;
    Redist = 0*D;
        Advect = 0*D;
    Swell = 0*D;
    DA = ones(size(D));
    t = 2;
    while (t < 25000) && (norm(DA) > tol)
        [A1,A2,V1,V2,epsdot] = load_bc(D,0);
        DA = 0*D;
    % Compute and updat the ice velocity and stress/strain tensors
    %% [Sigma,Epsdot,Vel,H] = do_ice_dynamics(t,H,Vel,Sigma,Epsdot,Voce,Vair);
    
   

    %% Advect Ice - Coded completely, no issues
      Advect = x_advect_fsd(A,D,V1,V2,A1,A2);       %.5 just because    
         Aloss = sum(Advect);
    
    
    
    %% Thermodynamic Melting - Coded
      Melt = melt_fsd(A,T,D);
    %Melt = 0;
    
    
    
    
    %% Mechanical Redistribution - Coded completely, no issues
    Redist = redist_fsd(A,epsdot,f,D,shiftra,shiftri);
    %[Part Raft Ridge] = Redist_OW_fsd(A,epsdot,f,D,shiftra,shiftri);
    Redist(1) = (Redist(1) - Aloss); % This handles the divergence
    
    %Redist = Part + Raft + Ridge; 
    
    
    %% Swell Fracture - Coded, Runs, maybe not ready
    
    Swell = swellfrac_fsd(A,D,lambda,Yg,target,Spectrum);

    

    %% We're Done!
    DA = Advect + Melt + Redist + Swell;
    Anew = A + dt*DA; 

    %% Compute relative magnitude of terms, stats, and bug control
    
    if min(Anew) < 0 
        disp(['Negative value at t = ' sprintf('%d',dt*t)])
        Ttotal = t;
        break
    else
        if isnan(max(Anew))
            disp(['NaN value at t = ' sprintf('%d',dt*t)])
            break
        else
            A = Anew;
        end
    end
    
    
 %    Save = Save + dt*[Advect; Melt; Redist; Swell];
    
 %   Save2 = Save2 + dt*Swell; 
    
%     if mod(t,dumpfreq) == 0
%         stats(t/dumpfreq,:) = [norm(Advect) norm(Melt) norm(Melt)/norm(Advect) ...
%          norm(Redist) norm(Redist)/norm(Advect) norm(Swell) ...
%          norm(Swell)/norm(Advect)];
%      
%         sumstats(t/dumpfreq,:) = [A(1) sum(Advect) sum(Melt) ...
%          sum(Redist) sum(Swell) sum(DA)];
%      sumstats2(t/dumpfreq,:) = Swell';
%      
%     Vals(t/dumpfreq,:,:) = [Melt; Redist; Advect;Swell];
    t = t + 1;
    end  
    
    Openwater(i) = A(1); 
    MeanD(i) = mean(A.*D);
    MeanA(i) = mean(A.*D.*D); 
    Max(i) = max(A);
    Dev(i) = std(A); 
end
% if norm(DA) < tol
%     disp(sprintf('Reached Steady State by t = %d',t*dt));
% end
stats = stats(1:t-1,:);
Vals = Vals(1:t-1,:,:);

%% Plot everything

saveplots = 0;

Saver = Save./(repmat(A0,[4 1])) + .0001;

Timer = linspace(0,dt*t,length(stats(1:t-1,1)));

plotFSD(D,[A;A0;A1],Timer', ... 
    [0*Timer + 1;stats(:,3)';stats(:,5)';stats(:,7)'], ... 
    [stats(:,1)';stats(:,2)';stats(:,4)';stats(:,6)'], ... 
    abs(Saver)/(dt*nt) + .0000001,squeeze(Vals(end,:,:)),t,saveplots);

%Tquit = round((t-1)/dumpfreq);

%plotbalance(D,[A;A0;A1],squeeze(Vals(end,:,:)),t);