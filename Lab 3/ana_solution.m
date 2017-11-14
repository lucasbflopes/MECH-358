function [ sum ] = ana_solution( x,t )

% This functions evaluates the analytical solution of the heat equation
% with precribed boundary conditions obtained in P1. Since the solution
% was written as a Fourier Series, it has been computed up to the 100th non
% zero term.

sum = 0;

for k = 1:100
    
    sum = sum + 8/(pi^3*(2*k-1)^3)*sin(pi*(2*k-1).*x).*...
        exp(-pi^2*(2*k-1)^2.*t);
    
end

end

