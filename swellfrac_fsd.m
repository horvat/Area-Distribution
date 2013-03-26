function DA = swellfrac_fsd(A,D,lambda,Yg,target,Spectrum)
% function DA = swellfrac_fsd(A,lambda,H)
% This function computes the change in the area distribution due to swell
% fracture, by incorporating wave attenuation effects and effects due to
% the thickness distribution

% A is the given area distribution
% D is a vector containing the leftmost area in each area box
% Tp is the dominant period for the Pierson-Moskovitz-Bretschneider Spectrum
% Hs is the significant wave height for the PMB spectrum
% H is the thickness distribution thicknesses
% g is the thickness distribution


%% Calculate the swell fracture part

%Calculate for each lambda, and each area category
Sp = Spectrum/sum(Spectrum);
numbins = length(D);
numl = length(lambda);
Loss = zeros(numl,numbins);
Gain = zeros(1,numbins);

%% 
Bmat = bsxfun(@gt,D',lambda);

for i = 1:numl
    
   Bvec = Bmat(:,i);
    
   Loss(i,:) = -(A.*Bvec')*(Yg(i)*Sp(i));

   for j = 1:numbins
   
       Gain(target(i)) = Gain(target(i)) + Yg(i)*Sp(i)*A(j)*(Bvec(j));
    
   end
        
    %Gain(i,:) = Yg(i)*Sp(i)*swellmode(A,D,target(i,:));
    
end

Loss = sum(Loss,1);
DA = Loss + Gain;


end

