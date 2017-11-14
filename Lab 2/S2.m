function [ Sn ] = S2( x,N )

% This function evaluates the partial sum of the Fourier Series of
% f(x) = { -1 on  -2<=x<=0        up to the Nth term.
%        {  1 on   0<=x<=2  

sum = 0;

for k = 1:N
    
    sum = sum + 4/((2*k-1)*pi)*sin((2*k-1)*pi*x/2);
end

Sn = sum;   % a0/2 = 0

end