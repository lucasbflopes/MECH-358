% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

function [ up ] = item_a(y,u)
% This function represents the RHS of the system of ODEs

G = -2;

x1 = u(1);
x2 = u(2);

up1 = x2;
up2 = G;

up=[up1 ; up2];

end
