function [ sum ] = ana_solution_P2( x,y )

% This functions evaluates the analytical solution of laplace's equation
% with precribed boundary conditions obtained in P2. Since the solution
% was written as a Fourier Series, it has been computed up to the 20th non
% zero term.

sum = 0;

for k = 1:20
    
    sum = sum - 2/(k*pi*sinh(k*pi))*((-1)^k - cos(k*pi/2))*...
        sinh(k*pi*x)*sin(k*pi*y);
    
end

end

