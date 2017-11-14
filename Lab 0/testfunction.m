% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

function [ f ] = testfunction( t,y )

% This function evaluates the RHS of dy/dt = f(t,y)

f = -3*y*t^2 + sin(t)*exp(-t^3);

end

