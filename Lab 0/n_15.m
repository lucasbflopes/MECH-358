% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

function [ up ] = n_15( t,y)
% This function represents the RHS of the system of ODEs

n = 15;

x1 = y(1);
x2 = y(2);
up1 = x2;
up2 = -x1^n;
up=[up1 ; up2];

end

