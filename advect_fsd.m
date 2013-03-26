function Advect = advect_fsd(A,Vel)
% This function computes the changes in the floe size distribution due to
% advection of the ice by the underlying flow

%A is the FSD being advected. 
%V is the 2D floe velocity field, defined on the Arakawa A-Grid for now.
%Vel(1:nx,1:ny,1) is the zonal velocity, Vel(1:nx,1:ny,2) is the 
%meridional velocity. 

%% Get Initial Conditions


% Gets boundary conditions for the problem. Defined on a grid enlarged by
% one. 

[Abc,Vel] = get_bc; 

% Gets grid size
[nx ny] = size(A); 

%Splits velocities
Uvel = Vel(:,:,1);
Vvel = Vel(:,:,2);

%% Set up adjusted A, V matrices

A = makenew(A,Abc);
Uvel = makenew(UVel,Ubc);
Vvel = makenew(VVel,Vbc);        



%% Advect in X
% Using upwind advection scheme (I think)
for i = 1:nx
   %Must account for extra padded size
   
   %figure out advection schemes, please. 

end

end


function Qnew = makenew(Q,Qbc)
% This function embeds the old matrix into a slightly larger one for flux
% computations
[nx,ny] = size(Q);
Qnew = zeros(nx+2,ny+2);
Qnew(1,:) = Qbc(1,:);
Qnew(nx+2,:) = Qbc(nx,:); 
Qnew(:,1) = Qbc(:,1);
Qnew(:,ny+2) = Qbc(:,ny+2);
Qnew(2:nx+1,2:ny+1) = Q; 
end