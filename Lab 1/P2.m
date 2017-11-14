% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve the BVP y''(x) = -[y'(x)]^2 - 0.1y + Ln(x) 
% y(1) = 0 , y(2) = 2
% In order to do so, it has been implemented a shooting method combined
% with bisection method.

% All instructions about what has been asked to do in this problem 
% can be found in Lab1.2015.pdf available at http://www.math.ubc.ca/
% ~coombs/358/math358.html

clear all ; clc

% RHS of the system   x1 = x2    , x1 = y(x)  and x2 = y'(x) 
%                     x2 = -x2^2 - 0.1x1 + Ln(x)  

f = @(x,y) [y(2);-y(2)^2-0.1*y(1)+log(x)];

% Initial guesses

a1 = -8;    % first guess of y'(1)
b1 =12;     % second guess of y'(1) 

% Calculates y(2) by ode45

[~,a_1] = ode45(f,[1 2],[0,a1]);
a11 = a_1(length(a_1(:,1)),1);   % y(2) when y'(1) = first guess
[~,b_1] = ode45(f,[1 2],[0,b1]);
b11 = b_1(length(b_1(:,1)),1);   % y(2) when y'(1) = second guess

% Creates a function g to tell how much close we are from the
% solution y(2) = 2

g = @(x) x - 2;    

% Starts the Shooting method combined with the Bisection method

if g(a11)*g(b11) > 0
    display('Your guesses are not valid. Try other ones')
    return
end   

c1 = (a1+b1)/2 ;  % estimated y'(1)
[x,u_num] = ode45(f,[1 2],[0,c1]);
c11 = u_num(length(u_num(:,1)),1);   % estimated y(2)
iter = 1;   % # of iterations

derivative_guess(iter) = c1;
finalvalue_guess(iter) = c11;
relative_error(iter) = [(c11-2)/2]*100;

e = 0.001; % absolute error ( precision up to the 3th decimal place)
while abs(g(c11)) > e 

    if g(c11)*g(a11) < 0
        b1 = c1;
        b11 = c11;
    else
        a1 = c1;
        a11 = c11;
    end

    c1 = (a1+b1)/2 ; 
    [x,u_num] = ode45(f,[1 2],[0,c1]);
    c11 = u_num(length(u_num(:,1)),1);

    iter = iter+1;

    derivative_guess(iter) = c1;
    finalvalue_guess(iter) = c11;
    relative_error(iter) = [(c11-2)/2]*100;

end

% Shows y'(1) obtained numerically

display(['The numerical value obtained for y ''(1) = ',num2str(c1)])

% Plot required curves

figure(1); 
plot(x,u_num(:,1)) 
xlabel('x');ylabel('y(x)')
title('Numerical solution for the BVP')

figure(2); 
plot(derivative_guess,'r-')
hold on
plot(finalvalue_guess,'b-')
xlabel('# of iterations');
set(gca,'XTick',1:1:length(derivative_guess))
title('Plot of estimates vs iteration')
legend('y ''(1)','y(2)');
hold off

figure(3);
plot(relative_error)
xlabel('# of iterations'); ylabel('Relative error (%)');
title('Relative error in y(2) vs # of iterations');


