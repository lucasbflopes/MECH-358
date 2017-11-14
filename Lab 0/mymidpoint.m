% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

% This file implements Midpoint method for the ODE in Problem 2

% External m.files used in this code:
% -> testfunction.m

clear all 
clc

% Initial conditions

t0 = 0;
y0 = 1;

% Other useful information

tfinal = 1; % is given
h = 0.1; % step size
N = (tfinal-t0)/h;  % number of points between t0 and tfinal

y = zeros(N+1,1);
t = zeros(N+1,1);
y(1) = y0;
t(1) = t0;

for i = 1:N
    y(i+1) = y(i) + h*testfunction(t(i)+h/2,y(i)+h/2*testfunction(...
        t(i),y(i)));
    t(i+1) = t(1) + i*h;
end


    