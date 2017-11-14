% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

% This file uses Matlab ode45 function to solve the ODE in Problem 2

% External m.files used in this code:
% -> testfunction.m

clear all
clc

y0 = 1;   % initial condition

% Call function

[t,y] = ode45(@testfunction,[0 1],y0);


