function [ x,D2_f ] = cen2( f,a,b,h )

% This function calculates the second derivative of the function f by the 
% central difference method on a range a <= x <= b, with step size h

%------ Central difference scheme------

%   f''(x) = [f(x+h) - 2f(x) + f(x-h)]/h^2 + O(h^2)

N = (b-a)/h ;

D2_f = zeros(1,N-1);
x = zeros(1,N-1);

for i = 1:N-1
    
    x(i) = a + i*h;
    D2_f(i) = (f(x(i)+h) - 2*f(x(i)) + f(x(i)-h))/h^2;

end

end

