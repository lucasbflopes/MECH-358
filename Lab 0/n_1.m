% Name: Lucas Bragan�a Fernandes lopes
% Student number: 56412142

function [ up ] = n_1( t,y)
% This function represents the RHS of the system of ODEs

n = 1;

x1 = y(1);
x2 = y(2);
up1 = x2;
up2 = -x1^n;
up=[up1 ; up2];

end

