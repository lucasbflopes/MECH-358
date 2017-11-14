function [ x,D_f ] = fwd1( f,a,b,h )

% This function calculates the derivative of the function f by the forward
% difference method on an interval a <= x <= b, with step size h

%------ Forward difference scheme------

%   f'(x) = [f(x+h)-f(x)]/h + O(h)

N = (b-a)/h ;

D_f = zeros(1,N);
x = zeros(1,N);

for i = 1:N
    
    x(i) = a + (i-1)*h;
    D_f(i) = (f(x(i)+h)-f(x(i)))/h;

end

end

