function [ x,D_f ] = bwd1( f,a,b,h )

% This function calculates the derivative of the function f by the backward
% difference method on an interval a <= x <= b, with step size h

%------ Backward difference scheme------

%   f'(x) = [f(x)-f(x-h)]/h + O(h)

N = (b-a)/h ;

D_f = zeros(1,N);
x = zeros(1,N);

for i = 1:N
    
    x(i) = a + i*h;
    D_f(i) = (f(x(i))-f(x(i)-h))/h;

end

end

