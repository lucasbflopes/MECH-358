function [ sum ] = u1( r,z,zeros)

% Analytical solution of the problem with the boundary condition:
% u(r,0) = 100

sum = 0;

for k = 1:length(zeros)
    
    sum = sum -200/(tanh(zeros(k))*zeros(k)*besselj(1,zeros(k)))*...
        besselj(0,zeros(k)*r).*(sinh(zeros(k)*z) - tanh(zeros(k))*...
        cosh(zeros(k)*z));

end

