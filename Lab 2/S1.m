function [ Sn ] = S1( x,N )

% This function evaluates the partial sum of the Fourier Series of
% f(x) = |x| on  -2<=x<=2 up to the Nth term.

sum = 0;

for k = 1:N
    
    sum = sum -8/((2*k-1)^2*pi^2)*cos((2*k-1)*pi*x/2);
    
end

Sn = 1 + sum;   % a0/2 = 1

end

