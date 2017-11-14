% Name: Lucas Bragança Fernandes lopes
% Student number: 56412142

% This code attempts to solve by using three different methods the follow-
% ing ODE:  y'=-3yt^2 + sin(t)*exp(-t^3), y(0) = 0 , 0 <= t <= 1
% These are: Euler's method, Midpoint method and Matlab function ode45
% Furthermore, all the steps showed on the instructions are followed here

% External m.files used in this code:
% -> myeuler.m
% -> mymidpoint.m
% -> myode45solver.m
% -> testfunction.m

clear all
clc

% This first part is responsible for calculating the solution for each
% method and then plotting the results in a single graph.

t0 = 0;
y0 = 1;

tfinal = 1; % is given
h = 0.1; % step size
N = (tfinal-t0)/h;  % number of points between t0 and tfinal

y1 = zeros(N+1,1);
t1 = zeros(N+1,1);
y1(1) = y0;
t1(1) = t0;

y2 = zeros(N+1,1);
t2 = zeros(N+1,1);
y2(1) = y0;
t2(1) = t0;

% Euler's method

for i = 1:N
    y1(i+1) = y1(i) + h*testfunction(t1(i),y1(i));
    t1(i+1) = t1(1) + i*h;
end

% Midpoint method

for i = 1:N
    y2(i+1) = y2(i) + h*testfunction(t2(i)+h/2,y2(i)+h/2*testfunction(...
        t2(i),y2(i)));
    t2(i+1) = t2(1) + i*h;
end

% ode45 method

[t3,y3] = ode45(@testfunction,[0 1],y0);

% One can easily show that the analytical solution for the ODE is:

y_analitical_1 = exp(-t1.^3).*(2-cos(t1));   
y_analitical_2 = exp(-t2.^3).*(2-cos(t2));
y_analitical_3 = exp(-t3.^3).*(2-cos(t3));

% The absolute error of the methods are:

error1 = abs(y1 - y_analitical_1);
error2 = abs(y2 - y_analitical_2);
error3 = abs(y3 - y_analitical_3);

% Ploting the erros in a single graph

figure;
semilogy(t1,error1,'r-+')
hold on
semilogy(t2,error2,'g-v')
semilogy(t3,error3,'b-s')
xlabel('t')
ylabel('absolute error')
legend('Euler','Midpoint','ode45','Location','northwest')
title('Absolute error in each method')
hold off

% The second part is responsible for calculating the solution of the ODE
% through Euler' method for a range a different step sizes (h)

clear all
clc

h = [0.2 0.1 0.05 0.025 0.0125];

t0 = 0;
y0 = 1;
tfinal = 1;

% h = 0.2
 
N1 = (tfinal-t0)/h(1);
y1 = zeros(N1+1,1);  
t1 = zeros(N1+1,1);    
y1(1) = y0;
t1(1) = t0;
for i = 1:N1
    y1(i+1) = y1(i) + h(1)*testfunction(t1(i),y1(i));
    t1(i+1) = t1(1) + i*h(1);
end

y_analitical_1 = exp(-t1.^3).*(2-cos(t1));

error1 = abs(y1 - y_analitical_1);

% h = 0.1

N2 = (tfinal-t0)/h(2);
y2 = zeros(N2+1,1);  
t2 = zeros(N2+1,1);    
y2(1) = y0;
t2(1) = t0;
for i = 1:N2
    y2(i+1) = y2(i) + h(2)*testfunction(t2(i),y2(i));
    t2(i+1) = t2(1) + i*h(2);
end

y_analitical_2 = exp(-t2.^3).*(2-cos(t2));

error2 = abs(y2 - y_analitical_2);

% h = 0.05

N3 = (tfinal-t0)/h(3);
y3 = zeros(N3+1,1);  
t3 = zeros(N3+1,1);    
y3(1) = y0;
t3(1) = t0;
for i = 1:N3
    y3(i+1) = y3(i) + h(3)*testfunction(t3(i),y3(i));
    t3(i+1) = t3(1) + i*h(3);
end

y_analitical_3 = exp(-t3.^3).*(2-cos(t3));

error3 = abs(y3 - y_analitical_3);

% h = 0.025
 
N4 = (tfinal-t0)/h(4);
y4 = zeros(N4+1,1);  
t4 = zeros(N4+1,1);    
y4(1) = y0;
t4(1) = t0;
for i = 1:N4
    y4(i+1) = y4(i) + h(4)*testfunction(t4(i),y4(i));
    t4(i+1) = t4(1) + i*h(4);
end

y_analitical_4 = exp(-t4.^3).*(2-cos(t4));

error4 = abs(y4 - y_analitical_4);

% h = 0.0125

N5 = (tfinal-t0)/h(5);
y5 = zeros(N5+1,1);  
t5 = zeros(N5+1,1);    
y5(1) = y0;
t5(1) = t0;
for i = 1:N5
    y5(i+1) = y5(i) + h(5)*testfunction(t5(i),y5(i));
    t5(i+1) = t5(1) + i*h(5);
end

y_analitical_5 = exp(-t5.^3).*(2-cos(t5));

error5 = abs(y5 - y_analitical_5);

figure;
semilogy(t1,error1,'r-*')
hold on                 % plot everything
semilogy(t2,error2,'g-o')
semilogy(t3,error3,'k-+')
semilogy(t4,error4,'b-v')
semilogy(t5,error5,'m-s')
xlabel('t')
ylabel('absolute error')
title('Absolute error in Eulers method for different step sizes')
legend('h = 0.2','h = 0.1','h = 0.05','h = 0.025','h = 0.0125',...
    'Location','northwest')
hold off

% Error at t=1 as a function of h

error_t_1 =[error1(length(error1)),error2(length(error2))...
    error3(length(error3)),error4(length(error4)),error5(length(error5))];

figure;
plot(h,error_t_1)
xlabel('h')
ylabel('absolute error')
title('Absolute error at t = 1 as function of step size (h) [Eulers method]')


% The third part does the same thing as second part but for the Midpoint
% method.

clear all
clc

h = [0.2 0.1 0.05 0.025 0.0125];

t0 = 0;
y0 = 1;
tfinal = 1;

% h = 0.2
 
N11 = (tfinal-t0)/h(1);
y11 = zeros(N11+1,1);  
t11 = zeros(N11+1,1);    
y11(1) = y0;
t11(1) = t0;
for i = 1:N11
    y11(i+1) = y11(i)+h(1)*testfunction(t11(i)+h(1)/2,y11(i)+h(1)/2*...
        testfunction(t11(i),y11(i)));
    t11(i+1) = t11(1) + i*h(1);
end

y_analitical_11 = exp(-t11.^3).*(2-cos(t11));

error11 = abs(y11 - y_analitical_11);

% h = 0.1

N22 = (tfinal-t0)/h(2);
y22 = zeros(N22+1,1);  
t22 = zeros(N22+1,1);    
y22(1) = y0;
t22(1) = t0;
for i = 1:N22
    y22(i+1) = y22(i)+ h(2)*testfunction(t22(i)+h(2)/2,y22(i)+h(2)/2*...
        testfunction(t22(i),y22(i)));
    t22(i+1) = t22(1) + i*h(2);
end

y_analitical_22 = exp(-t22.^3).*(2-cos(t22));

error22 = abs(y22 - y_analitical_22);

% h = 0.05

N33 = (tfinal-t0)/h(3);
y33 = zeros(N33+1,1);  
t33 = zeros(N33+1,1);    
y33(1) = y0;
t33(1) = t0;
for i = 1:N33
    y33(i+1) = y33(i)+h(3)*testfunction(t33(i)+h(3)/2,y33(i)+h(3)/2*...
        testfunction(t33(i),y33(i)));
    t33(i+1) = t33(1) + i*h(3);
end

y_analitical_33 = exp(-t33.^3).*(2-cos(t33));

error33 = abs(y33 - y_analitical_33);

% h = 0.025
 
N44 = (tfinal-t0)/h(4);
y44 = zeros(N44+1,1);  
t44 = zeros(N44+1,1);    
y44(1) = y0;
t44(1) = t0;
for i = 1:N44
    y44(i+1) = y44(i)+h(4)*testfunction(t44(i)+h(4)/2,y44(i)+h(4)/2*...
        testfunction(t44(i),y44(i)));
    t44(i+1) = t44(1) + i*h(4);
end

y_analitical_44 = exp(-t44.^3).*(2-cos(t44));

error44 = abs(y44 - y_analitical_44);

% h = 0.0125

N55 = (tfinal-t0)/h(5);
y55 = zeros(N55+1,1);  
t55 = zeros(N55+1,1);    
y55(1) = y0;
t55(1) = t0;
for i = 1:N55
    y55(i+1) = y55(i)+h(5)*testfunction(t55(i)+h(5)/2,y55(i)+h(5)/2*...
        testfunction(t55(i),y55(i)));
    t55(i+1) = t55(1) + i*h(5);
end

y_analitical_55 = exp(-t55.^3).*(2-cos(t55));

error55 = abs(y55 - y_analitical_55);

figure;
semilogy(t11,error11,'r-*')
hold on                 % plot everything
semilogy(t22,error22,'g-o')
semilogy(t33,error33,'k-+')
semilogy(t44,error44,'b-v')
semilogy(t55,error55,'m-s')
xlabel('t')
ylabel('absolute error')
title('Absolute error in Midpoint method for different step sizes')
legend('h = 0.2','h = 0.1','h = 0.05','h = 0.025','h = 0.0125',...
    'Location','northwest')
hold off

% Last part: error at t=1 as a function of h

error_t_11 =[error11(length(error11)),error22(length(error22))...
    error33(length(error33)),error44(length(error44)),error55(length...
    (error55))];

figure;
plot(h,error_t_11)
xlabel('h')
ylabel('absolute error')
title('Absolute error at t = 1 as function of step size (h) [Midpoint method]')


    
