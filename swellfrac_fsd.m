function DA = swellfrac_fsd(A,D,lambda,Yg,target,Spectrum,dt)
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
       if target(i) == 1
       Gain(target(i)) = Gain(target(i)) + .75*Yg(i)*Sp(i)*A(j)*(Bvec(j));
       Gain(target(i) + 1) = Gain(target(i) + 1) + .25*Yg(i)*Sp(i)*A(j)*(Bvec(j));
       
       elseif target(i) == numbins
                 Gain(target(i) - 1) = Gain(target(i) - 1) + .25*Yg(i)*Sp(i)*A(j)*(Bvec(j));
                 Gain(target(i)) = Gain(target(i)) + .75*Yg(i)*Sp(i)*A(j)*(Bvec(j));
           
           else
                 Gain(target(i) - 1) = Gain(target(i) - 1) + .25*Yg(i)*Sp(i)*A(j)*(Bvec(j));
                 Gain(target(i)) = Gain(target(i)) + .5*Yg(i)*Sp(i)*A(j)*(Bvec(j));
                 Gain(target(i) + 1) = Gain(target(i) + 1) + .25*Yg(i)*Sp(i)*A(j)*(Bvec(j));
              
       end
    end
       
end


        
    %Gain(i,:) = Yg(i)*Sp(i)*swellmode(A,D,target(i,:));

Loss = sum(Loss,1);
DA = Loss + Gain;
% 
% 
% if max(DA(2:numbins)) ~= 0
% s1 = sum(DA(2:numbins));
% DA(2:numbins) = real(exp(smooth(log(DA(2:numbins)),3)'));
% s2 = sum(DA(2:numbins));
% DA(2:numbins) = (s1/s2)*DA(2:numbins);
% end

% for i = 2:numbins-1
%     qDA = .25*DA(i); 
%     if dt*(DA(i-1) + .25*DA(i)) + A(i-1) > 0
%     DA(i-1) = DA(i-1) + .25*DA(i); 
%     DA(i) = DA(i) - qDA;
%     end
%     if dt*(DA(i+1) + .25*DA(i)) + A(i) > 0
%     DA(i+1) = DA(i+1) + .25*DA(i);
%     
%     DA(i) = DA(i) - qDA;
%     end
% end
end

