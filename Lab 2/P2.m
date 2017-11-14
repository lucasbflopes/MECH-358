% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve a two point boundary value problem using finite
% difference and to compare the numerical results against an exact solution

% External m.files used in this code:
% -> finite_difference.m

clear all ; clc

% Definition of relevant parameters and boundary conditions

h = .1 ; % Define mesh spacing
Pe = 1 ; % Define Peclét Number
L = 10 ; % Define length of the plate

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

T = [T0 Tm' T_inf]';   % Complete solution including the values of T at 
% the boundary points

% Analytical solution of the problem

x = 0:h:L;
T_ana = exp(-Pe*x);

clear all; clc;

% An auxiliary function has been created to facilitate the process of
% comparing the analytical solution and the numerical one for a range of
% different values of h and L. This funtion basically follows the same 
% steps presented above

% For the case where we are varying L, h was kept equal to .1

[x_10,T_num_10,T_ana_10] = finite_difference(.1,10);
[x_50,T_num_50,T_ana_50] = finite_difference(.1,50);
[x_100,T_num_100,T_ana_100] = finite_difference(.1,100);

% Plot numerical solution vs analytical solution for each L

% It has been chosen a semilog plot to clarify the differences between both
% solutions.

figure(1);
semilogy(x_10,T_num_10,'r-',x_10,T_ana_10,'b-')
xlabel('x')
ylabel('T(x)')
title('Temperture vs distance for L = 10 and h = 0.1')
legend('Numerical solution','Analytical solution')
print('P2_1', '-dpng', '-r300');

figure(2);
semilogy(x_50,T_num_50,'r-',x_50,T_ana_50,'b-')
xlabel('x')
ylabel('T(x)')
title('Temperture vs distance for L = 50 and h = 0.1')
legend('Numerical solution','Analytical solution')
print('P2_2', '-dpng', '-r300');

figure(3);
semilogy(x_100,T_num_100,'r-',x_100,T_ana_100,'b-')
xlabel('x')
ylabel('T(x)')
title('Temperture vs distance for L = 100 and h = 0.1')
legend('Numerical solution','Analytical solution')
print('P2_3', '-dpng', '-r300');

clear all; clc

% For the case where we are varying h and L was kept equal to 100

[x_1,T_num_1,T_ana_1] = finite_difference(1,100);
[x_05,T_num_05,T_ana_05] = finite_difference(.5,100);
[x_01,T_num_01,T_ana_01] = finite_difference(.1,100);


% Plot numerical solution vs analytical solution for each L

% OBS: It has been chosen a semilog plot to clarify the difference between 
% both solutions.

figure(4);
semilogy(x_1,T_num_1,'r-'); hold on;
semilogy(x_05,T_num_05,'b-');
semilogy(x_01,T_num_01,'g-');
xlabel('x')
ylabel('T(x)')
legend('h = 1','h = 0.5','h = 0.1')
title('Temperture vs distance for different values of h and L = 100')
print('P2_4', '-dpng', '-r300');




