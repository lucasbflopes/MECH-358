% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to compute with Bessel Series

% External m.files used in this code:
% -> u1.m
% -> u2.m
% -> zerobess.m

% URL of the webpage on where the code that finds the zeroes of the bessel
% function were downloaded:

% http://www.mathworks.com/matlabcentral/fileexchange/26639-zerobess


clear all; clc;


% Since the Bessel series representation of the solution of the problem has
% 20 terms, we have to calculate the first 20 zeros of the first kind bessel
% function of order zero

% According to a paper found on the web, which will be attached to the
% present work and its URL referenced on the head of this very document,
% the first 20 zeros are:

bessel_zeros = zerobess('J',0,20);

N = length(bessel_zeros); 

h = 0.01; % step size
r = 0:h:1; 
z = 0:h:1;
M = length(r);

% The solution of the steady-state heat equation with the prescribed
% boundary conditions and such geometry is:
%
%
%
% u(r,z) = sum_(k=1:20) Ak*(J0(lambdak*r)*(sinh(lambdak*z) - 
%                                           tanh(lambdak)*cosh(lambdak*z))
%
%
% where Ak = -2/(tanh(lambdak)*J1^2(lambdak)*integral_0:1
%                                               f(r)*J0(lambdak*r)*r*dr
%
%

% ------------------------------------------------------------------------
% ------------------------  Part 1: f(r) = 100  --------------------------
% ------------------------------------------------------------------------

% For this simple boundary condition we are able to calculate the Ak's
% analytically without any numerical integration. 

% Graphs

figure(1); 
plot(r,u1(r,0.5,bessel_zeros),'LineWidth',1.2)
set(gca,'FontSize',12);
grid on;
xlabel('r')
ylabel('T(r,0.5)')
title('Analytical solution for boundary condition u(r,0) = 100')
print('P2_1', '-dpng', '-r300');


% ------------------------------------------------------------------------
% ---------------------  Part 2: f(r) = cos(pi*r/2)  ---------------------
% ------------------------------------------------------------------------

% Defining the boundary condition at u(r,z = 0)

f = @(r) cos(pi*r/2);  

% In order to calculate Ak, we have to compute an integral. The method that
% has been chosen to do so was the Simpson's rule

Ak = zeros(1,N);

% Starting Simpson's rule

for i = 1:N
    
    integrand = f(r).*besselj(0,bessel_zeros(i)*r).*r;  
    integral = h/3*( integrand(1) + 4*sum(integrand(2:2:M-1))...
        + 2*sum(integrand(3:2:M-2)) + integrand(M) ) ;
    Ak(i) = -2/(tanh(bessel_zeros(i))*besselj(1,bessel_zeros(i))^2)*...
        integral;
    
end

% Graphs

figure(2);
plot(r,u2(r,0.5,Ak,bessel_zeros),'LineWidth',1.2)
set(gca,'FontSize',12);
grid on;
xlabel('r')
ylabel('T(r,0.5)')
title('Numerical solution for boundary condition u(r,0) = cos(pi*r/2)')
print('P2_2', '-dpng', '-r300');

    
       
    
