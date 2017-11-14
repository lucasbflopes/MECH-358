function [ x,D_f ] = cen1( f,a,b,h )

% This function calculates the derivative of the function f by the central
% difference method on a range a <= x <= b, with step size h

%------ Central difference scheme------

%   f'(x) = [f(x+h)-f(x-h)]/2h + O(h^2)

N = (b-a)/h ;

D_f = zeros(1,N-1);
x = zeros(1,N-1);

for i = 1:N-1
    
    x(i) = a + i*h;  
    D_f(i) = (f(x(i)+h)-f(x(i)-h))/(2*h);

end

end

