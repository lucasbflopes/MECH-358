% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve the BVP u''(y) = 0 , u(-b) = 0 , u(b) = V
% In order to do so, it has been implemented a shooting method combined
% with bisection method.

% All instructions about what has been asked to do in item a and b 
% can be found in Lab1.2015.pdf available at http://www.math.ubc.ca/
% ~coombs/358/math358.html


%------------------------ Item a ----------------------------

clear all ; clc

b = 1;
V = 0;
G = -2;

% Initial condition

y_prime = 2;

% RHS of the system  x1 = x2      x1 = u(y), x1 = u'(y)
%                    x2 = G   

f = @(y,u) [u(2);G];

% Find solution by numerical method 
% (OBS: default relative error in ode45 is 10^-3)

[y,u_num] = ode45(f,[-b b],[0,y_prime]);  

% Solution obtained analytically

u_ana = 1 - y.^2;

% Plot both solutions to make a comparison

figure(1);
plot(y,u_ana,'r*',y,u_num(:,1),'b-')
xlabel('y'); ylabel('u(y)');
legend('Analytical solution','Numerical solution')
title('Plot of Numerical solution x Analytical solution for b=1,G=-2,V=0')

% ----------------------- Item b ----------------------------

clear all; clc;

b = 1;
V = 1;
G = -2;
e = 0.0001;   % absolute error ( precision up to the 4th decimal place)

% Initial guesses

a1 = 1;    % first guess of u'(-b) 
b1 = 4;    %  second guess of  u'(-b) 

% RHS of the system

f = @(y,u) [u(2);G];

% Calculate u(b) by ode45  

[~,a_1] = ode45(f,[-b b],[0,a1]);
a11 = a_1(length(a_1(:,1)),1);   % u(b) when u'(-b) = first guess
[~,b_1] = ode45(f,[-b b],[0,b1]);  
b11 = b_1(length(b_1(:,1)),1);   % u(b) when u'(-b) = second guess

% Creates a function g to tell how much close we are from the
% solution u(b) = V

g = @(x) x - V ;

% Starts the Shooting method combined with the Bisection method

if g(a11)*g(b11) > 0
    display('Your guesses are not valid. Try other ones')
    return
end   

c1 = (a1+b1)/2 ;  % estimated u'(-b)
[y,u_num] = ode45(f,[-b b],[0,c1]);
c11 = u_num(length(u_num(:,1)),1);  % estimated u(b)

while abs(g(c11)) > e 

    if g(c11)*g(a11) < 0
        b1 = c1;
        b11 = c11;
    else
        a1 = c1;
        a11 = c11;
    end

    c1 = (a1+b1)/2 ; 
    [y,u_num] = ode45(f,[-b b],[0,c1]);
    c11 = u_num(length(u_num(:,1)),1);

end

% Shows u'(-b) obtained numerically

display(['The numerical value obtained for u ''(-',num2str(b),') = '...
    ,num2str(c1)])

% Solution obtained analytically

u_ana = 1/2*(3 + y -2*y.^2);

% Plot both solutions to make a comparison

figure(2);
plot(y,u_ana,'r*',y,u_num(:,1),'b-')
xlabel('y'); ylabel('u(y)');
legend('Analytical solution','Numerical solution','Location','Northwest')
title('Plot of Numerical solution x Analytical solution for b=1,G=-2,V=1')





