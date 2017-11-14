function [ sum ] = ana_solution_P1( x,t )

% This functions evaluates the analytical solution of the heat equation
% with precribed boundary conditions obtained in P1. Since the solution
% was written as a Fourier Series, it has been computed up to the 100th non
% zero term.

sum = 0;

for k = 1:100  % 100 time steps
    
    sum = sum + 4/((2*k-1)*pi)*sin(pi*(2*k-1).*x).*...
        exp(-pi^2*(2*k-1)^2.*t);
    
end

end

