% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

% This code attempts to solve by using the Matlab function ode45 the second 
% order differential equation u''(t)+(u(t))^n = 0, u(0) = 0 and u'(0) = 1
% for n = 1,5,9,15

% External m.files used in this code:
% -> n_1.m
% -> n_5.m
% -> n_9.m
% -> n_15.m

clear all
clc

% Solve the associated ODE system for n = 1,5,9,15

[t_1,y_1] = ode45(@n_1,[0 20],[0;1]); % n = 1
[t_5,y_5] = ode45(@n_5,[0 20],[0;1]); % n = 5
[t_9,y_9] = ode45(@n_9,[0 20],[0;1]); % n = 9   
[t_15,y_15] = ode45(@n_15,[0 20],[0;1]);  % n = 15
y_1_analitical = sin(t_1);   % analitical solution when n = 1

% Compare both analytical and numerical solution

figure;
plot(t_1,y_1(:,1),'--ro')
hold on
plot(t_1,y_1_analitical,'--b*')  
xlabel('t')
ylabel('u(t)')
legend('Numerical solution','Analytical solution')
title('Comparison between analytical and numerical solution')
hold off

% Plot the solution for all n's in one plot

figure;
plot(t_1,y_1(:,1),'r') 
hold on
plot(t_5,y_5(:,1),'g')
plot(t_9,y_9(:,1),'b')
plot(t_15,y_15(:,1),'k')
xlabel('t')
ylabel('u(t)')
legend('n=1','n=5','n=9','n=15')
title('Solution of the ODE for different values of n')
hold off

% Energy equation for n = 1,5,9,15 :

energy_1 = (1/2)*y_1(:,2).^2 + (1/2)*y_1(:,1).^2; % n = 1
energy_5 = (1/2)*y_5(:,2).^2 + (1/6)*y_5(:,1).^6; % n = 5
energy_9 = (1/2)*y_9(:,2).^2 + (1/10)*y_9(:,1).^10; % n = 9
energy_15 = (1/2)*y_15(:,2).^2 + (1/16)*y_15(:,1).^16; % n = 15

% Since we know the analytical solution for n = 1, if we plug it in the 
% energy equation one might find out that C = 0.5. So, as long as C
% remains ~ 0.5, the values obtained numericaly are fine.

figure;
plot(t_1,energy_1,'r')
hold on
plot(t_5,energy_5,'g')
plot(t_9,energy_9,'b')
plot(t_15,energy_15,'k')
xlabel('t')
ylabel('LHS of energy equation')
legend('n=1','n=5','n=9','n=15','Location','northwest')
title('LHS of energy equation for different values of n')
hold off




