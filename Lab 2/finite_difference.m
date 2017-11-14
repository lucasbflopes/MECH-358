function [ x,T_num,T_ana ] = finite_difference( h1,L1 )

% This function is nothing but the algorithm developed in P2.m . Its
% purpose is to find out the temperature distribution along the plate
% for different values of h and L 

% Definition of relevant parameters and boundary conditions

h = h1 ; % Define mesh spacing
Pe = 1 ; % Peclé number is assumed to be constant and equals to 1
L = L1 ; % Define length of the plate

N = L/h - 1 ; % # of points between the ends of the mesh

T0 = 1;      % T(x = 0) = 1
T_inf = 0;   % T(x = infinity) = 0 

% OBS: Since in the numerical problem the plate is finite we use the
% boundary condition T(x = L) = 0 instead of the aforementioned one.

% Building A matrix ( coefficient matrix )

diagonal = 2*ones(1,N);

superdiagonal = (-1-h/2*Pe)*ones(1,N-1);

subdiagonal = (-1+h/2*Pe)*ones(1,N-1);

A = diag(diagonal) + diag(superdiagonal,1) + diag(subdiagonal,-1); 

% Building b vector ( right hand side vector )

b = zeros(N,1);
b(1) = (1-h/2*Pe)*T0;
b(N) = (1+h/2*Pe)*T_inf;

% Solve for T

Tm = A\b;   % Solution of the system

T_num = [T0 Tm' T_inf]';   % Complete solution including the values of T at 
% the boundary points

% Analytical solution of the problem

x = 0:h:L;
T_ana = exp(-Pe*x);

end

