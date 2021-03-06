function DA = redist_fsd_v2(A,epsdot,f,D,shiftra,shiftri,epsri,epsra,Multri,Multra)

%This function controls the binned redistribution, based on new ridging
%modes different from Thorndike and others

%The redistribution function looks like:
% Psi = |epsdot|*(A*delta(d) + B*(-P(r)a(r) + f*w_ra(r) +(1-f)*w_ri(r)))

% w_ra and w_ri are in fact integrals 

% w_ra = int_d1 int_d2 P(d1)P(d2)a(d1)a(d2)delta(d1+d2-eps_ra)

% Meaning we are looking for two partners whose new ridged/rafted area will
% match that of d. Before we just assumed a partner was "found". 

% Shiftrav2 and Shiftriv2 are matrices which describe which box d_i and d_j
% will add their area to in a given interaction

    % The assumption is: two pieces come together
    % The smallest of the two has some part of it turned into the ridge
    % (This is no different than saying both ridge, since they'll combine)
    % r percent of the original piece ridges/rafts
    % This turns into a ridge of height k times the original height
    % Therefore the original piece has an area which is now: 
    % A1 = (1-r)A1 + rA1 --> (1-r)A1 + (r/k)A1 = (1 - (k-1)r/k)A1 = A1'
    % And the total area becomes A2+A1 --> A2 + A1'
    % Plus an open water term
    
 % In ridging, k ~ 5 and r ~ .5 So half of the incident piece is
 % turned into a ridge
 % In rafting, k ~ 2 and r ~ 1 So the entire piece is turned into a double
 % thick piece. Or, the whole piece uplifts onto the other piece and half
 % of it remains. Same idea. 
 
    % This gives for the new diameter:
    % dnew = (d1^2 + (1 - (k-1)r/k)d2^2)^(1/2)

%% Basic Instantiation Stuff
mult = 100; %Artificial multiplier incase we want to boost redist
numbins = length(A);
inv = eig(epsdot);
% Strain Rate Invariants - Properly normalized
% Normalization doesn't affect the closing/opening ratio, thankfully
eps1 = (1/sqrt(2))*(inv(1) + inv(2)); %First strain rate invariant: the divergence
eps2 = (1/sqrt(2))*(inv(1) - inv(2)); %Second strain rate invariant: ~ the shear
strainmag = norm(epsdot,'fro'); %Frobenius Norm = root of sum of squares of entries
theta = atan(eps2/eps1); %Ratio of shear to divergence
%Opening and closing coefficients
leadopen = .5*(1 + cos(2*theta));
leadclose = .5*(1 - cos(2*theta));

%% Participation Kernel

% Use proposal.pdf to understand this equation


minn = bsxfun(@min,(A./(pi/4*D.^2))',A./(pi/4*D.^2)); %Number per sq meter

epsri = D/5;
epsra = D/2;

top = 16*(A)'*(A); 

% Top may have to include a A/D part 

ari = A.*(1 - epsri./D); 
ara = A.*(1 - epsra./D);

bottri = 1 - bsxfun(@minus,ari,ari');
bottra = 1 - bsxfun(@minus,ara,ara'); 

Particra = triu(min(top./bottra,minn));
Particri = triu(min(top./bottri,minn));

% 
% Partic = 16*(A.*A)
% Partic = triu(min(2*(A.*D)'*(A./(D.*(1 - A).^2)),minn));
 Particra(1,:) = 0; 
 Particra(:,1) = 0;
 Particri(:,1) = 0;
 Particri(1,:) = 0;
 
 C = max(max(max(Particra)),max(max(Particri))); 

if C > 1
    disp(C)
end
%% Redistribution

DA(1) = leadopen;

DAraft = 0*D;
DAridge = 0*D; 
DAr = DAridge;
%% Ridging and Rafting Mode
Partic = f*Particra + (1-f)*Particri; 

DAr = -sum(Partic,2)'; % Loss for each diameter

% Would like to find a better way to do this

for i = 2:numbins
    for j = i:numbins
        indra = shiftra(i,j);
        indri = shiftri(i,j);
                
        DAraft(indra) = DAraft(indra) + f*Particra(i,j)*Multra(i,j);
        DAraft(1) = DAraft(1) + f*Particra(i,j)*(1 - Multra(i,j));
        DAridge(indri) = DAridge(indri) + (1-f)*Particri(i,j)*Multri(i,j);
        DAridge(1) = DAridge(1) + (1-f)*Particri(i,j)*(1 - Multri(i,j));
        
        % Some of the incident area is lost to open water formation
        
    end
end


DApre = DAr + DAraft + DAridge; 

%%
%DApre = horvsmooth(DApre);

% The real DA
DA = mult*strainmag*leadclose*DApre;

%% If we have no open water, redistribution must be adjusted so it stays that way

% DA(1) = -sum(DA(2:numbins));

%%

% [C,I] = min(A + DA);
% 
% if C < 0 
%     DAo = DA(1);
%     DA(1) = DAo - C; 
%     frac = abs(DA(1)/DAo);
%     DA(2:numbins) = frac*DA(2:numbins);
% end

 end

