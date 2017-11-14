function [ sol ] = u2( r,z,Ak,bessel_zeros )

% Analytical solution of the problem with the boundary condition:
% u(r,0) = cos(pi*r/2)

sol = 0;

for i = 1:20
    
    sol = sol + Ak(i).*besselj(0,bessel_zeros(i)*r).*(sinh(bessel_zeros(i)*z) - ...
        tanh(bessel_zeros(i)).*cosh(bessel_zeros(i)*z)); 

end

