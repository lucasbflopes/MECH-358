% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to compare the accuracy of different kinds of numerical
% approximations to derivatives. For the first derivative it has been used
% the forward-, backward- and central-difference, and for the second
% derivative just the central-difference. Two sets of functions has been
% used as references to compute the accuracy of each method.

% External m.files used in this code:
% -> testfunction.m
% -> testfunction2.m
% -> fwd1.m
% -> bwd1.m
% -> cen1.m
% -> cen2.m

%----------------------------------------------------------------
%------------------------ First Part ----------------------------
%------------------ for f(x) = sin(sqrt(x)) ---------------------
%----------------------------------------------------------------

clear all ; clc

h = [ 1 .1 .01 .001];  % step size array

a = 0;   
b = 10;

% Derivatives obtained through each method when h = 1

[xf_1,yf_1] = fwd1(@testfunction,a,b,h(1));
[xb_1,yb_1] = bwd1(@testfunction,a,b,h(1));
[xc1_1,yc1_1] = cen1(@testfunction,a,b,h(1));
[xc2_1,yc2_1] = cen2(@testfunction,a,b,h(1));

% Derivatives obtained through each method when h = 0.1

[xf_2,yf_2] = fwd1(@testfunction,a,b,h(2));
[xb_2,yb_2] = bwd1(@testfunction,a,b,h(2));
[xc1_2,yc1_2] = cen1(@testfunction,a,b,h(2));
[xc2_2,yc2_2] = cen2(@testfunction,a,b,h(2));

% Derivatives obtained through each method when h = 0.01

[xf_3,yf_3] = fwd1(@testfunction,a,b,h(3));
[xb_3,yb_3] = bwd1(@testfunction,a,b,h(3));
[xc1_3,yc1_3] = cen1(@testfunction,a,b,h(3));
[xc2_3,yc2_3] = cen2(@testfunction,a,b,h(3));

% Derivatives obtained through each method when h = 0.001

[xf_4,yf_4] = fwd1(@testfunction,a,b,h(4));
[xb_4,yb_4] = bwd1(@testfunction,a,b,h(4));
[xc1_4,yc1_4] = cen1(@testfunction,a,b,h(4));
[xc2_4,yc2_4] = cen2(@testfunction,a,b,h(4));

% First derivative obtained analytically

x = a:.01:b;
y_prime = cos(sqrt(x))./(2*sqrt(x));

% Second derivative obtained analytically

y_doubleprime = -cos(sqrt(x))./(4*x.^(3/2)) - sin(sqrt(x))./(4*x);

% Absolute error for each function for each value of h

% h = 1

errorf_1 = abs(yf_1 - cos(sqrt(xf_1))./(2*sqrt(xf_1)));
errorb_1 = abs(yb_1 - cos(sqrt(xb_1))./(2*sqrt(xb_1)));
errorc1_1 = abs(yc1_1 - cos(sqrt(xc1_1))./(2*sqrt(xc1_1)));
errorc2_1 = abs(yc2_1 + cos(sqrt(xc2_1))./(4*xc2_1.^(3/2)) + sin(...
    sqrt(xc2_1))./(4*xc2_1));

% h = 0.1

errorf_2 = abs(yf_2 - cos(sqrt(xf_2))./(2*sqrt(xf_2)));
errorb_2 = abs(yb_2 - cos(sqrt(xb_2))./(2*sqrt(xb_2)));
errorc1_2 = abs(yc1_2 - cos(sqrt(xc1_2))./(2*sqrt(xc1_2)));
errorc2_2 = abs(yc2_2 + cos(sqrt(xc2_2))./(4*xc2_2.^(3/2)) + sin(...
    sqrt(xc2_2))./(4*xc2_2));

% h = 0.01

errorf_3 = abs(yf_3 - cos(sqrt(xf_3))./(2*sqrt(xf_3)));
errorb_3 = abs(yb_3 - cos(sqrt(xb_3))./(2*sqrt(xb_3)));
errorc1_3 = abs(yc1_3 - cos(sqrt(xc1_3))./(2*sqrt(xc1_3)));
errorc2_3 = abs(yc2_3 + cos(sqrt(xc2_3))./(4*xc2_3.^(3/2)) + sin(...
    sqrt(xc2_3))./(4*xc2_3));

% h = 0.001

errorf_4 = abs(yf_4 - cos(sqrt(xf_4))./(2*sqrt(xf_4)));
errorb_4 = abs(yb_4 - cos(sqrt(xb_4))./(2*sqrt(xb_4)));
errorc1_4 = abs(yc1_4 - cos(sqrt(xc1_4))./(2*sqrt(xc1_4)));
errorc2_4 = abs(yc2_4 + cos(sqrt(xc2_4))./(4*xc2_4.^(3/2)) + sin(...
    sqrt(xc2_4))./(4*xc2_4));

% Plot the numerical solution for each step size along with the analytical
% solution

% Forward-Difference

figure(1);
plot(x,y_prime,'r-'); hold on
plot(xf_1,yf_1,'b-')
plot(xf_2,yf_2,'g-')
plot(xf_3,yf_3,'c-')
plot(xf_4,yf_4,'k-')
ylim([-0.5 2])
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [FWD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_1', '-dpng', '-r300');

% Backward-Difference

figure(2);
plot(x,y_prime,'r-'); hold on
plot(xb_1,yb_1,'b-')
plot(xb_2,yb_2,'g-')
plot(xb_3,yb_3,'c-')
plot(xb_4,yb_4,'k-')
ylim([-0.5 2])
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [BWD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_2', '-dpng', '-r300');

% Central-Difference

figure(3);
plot(x,y_prime,'r-'); hold on
plot(xc1_1,yc1_1,'b-')
plot(xc1_2,yc1_2,'g-')
plot(xc1_3,yc1_3,'c-')
plot(xc1_4,yc1_4,'k-')
ylim([-0.5 2])
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [CD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_3', '-dpng', '-r300');

% Central-Difference (Second derivative)

figure(4);
plot(x,y_doubleprime,'r-'); hold on
plot(xc2_1,yc2_1,'b-')
plot(xc2_2,yc2_2,'g-')
plot(xc2_3,yc2_3,'c-')
plot(xc2_4,yc2_4,'k-')
ylim([-1.5 0.1])
xlabel('x')
ylabel('f''''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [CD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001','Location','Southeast')
print('P1_4', '-dpng', '-r300');

% Plot of the error in the derivative at x = 2 obtained from each 
% method vs h

% Forward-difference

figure(5);
loglog([h(1),h(2),h(3),h(4)],[errorf_1(3),errorf_2(21),errorf_3(201),...
    errorf_4(2001)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [FWD method]')
print('P1_5', '-dpng', '-r300');

% Backward-difference

figure(6);
loglog([h(1),h(2),h(3),h(4)],[errorb_1(2),errorb_2(20),errorb_3(200),...
    errorb_4(2000)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [BWD method]')
print('P1_6', '-dpng', '-r300');

% Central-difference

figure(7);
loglog([h(1),h(2),h(3),h(4)],[errorc1_1(2),errorc1_2(20),errorc1_3(200),...
    errorc1_4(2000)])
xlabel('h')
ylabel('absolute error ')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [CD method]')
print('P1_7', '-dpng', '-r300');

% Central-difference (Second Derivative)

figure(8);
loglog([h(1),h(2),h(3),h(4)],[errorc2_1(2),errorc2_2(20),errorc2_3(200),...
    errorc2_4(2000)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''''(x) at x = 2 vs step size (h)  [CD method]')
print('P1_8', '-dpng', '-r300');

% Plot of absolute error on the range

% Forward-difference

figure(9);
semilogy(xf_1,errorf_1,'r-',xf_2,errorf_2,'b-',xf_3,errorf_3,'g-',...
    xf_4,errorf_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [FWD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_9', '-dpng', '-r300');

% Backward-difference

figure(10);
semilogy(xb_1,errorb_1,'r-',xb_2,errorb_2,'b-',xb_3,errorb_3,'g-',...
    xb_4,errorb_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [BWD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_10', '-dpng', '-r300');

% Central-difference

figure(11);
semilogy(xc1_1,errorc1_1,'r-',xc1_2,errorc1_2,'b-',xc1_3,errorc1_3,'g-',...
    xc1_4,errorc1_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [CD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_11', '-dpng', '-r300');

% Central-difference (Second Derivative)

figure(12);
semilogy(xc2_1,errorc2_1,'r-',xc2_2,errorc2_2,'b-',xc2_3,errorc2_3,'g-',...
    xc2_4,errorc2_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''''(x) [CD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_12', '-dpng', '-r300');


%-----------------------------------------------------------------
%------------------------ Second Part ----------------------------
%------------------ for f(x) = x^2 - x + 1 -----------------------
%-----------------------------------------------------------------

clear all ; clc

h = [ 1 .1 .01 .001];  % step size array

a = -3;   
b = 3;

% Derivatives obtained through each method when h = 1

[xf_1,yf_1] = fwd1(@testfunction2,a,b,h(1));
[xb_1,yb_1] = bwd1(@testfunction2,a,b,h(1));
[xc1_1,yc1_1] = cen1(@testfunction2,a,b,h(1));
[xc2_1,yc2_1] = cen2(@testfunction2,a,b,h(1));

% Derivatives obtained through each method when h = 0.1

[xf_2,yf_2] = fwd1(@testfunction2,a,b,h(2));
[xb_2,yb_2] = bwd1(@testfunction2,a,b,h(2));
[xc1_2,yc1_2] = cen1(@testfunction2,a,b,h(2));
[xc2_2,yc2_2] = cen2(@testfunction2,a,b,h(2));

% Derivatives obtained through each method when h = 0.01

[xf_3,yf_3] = fwd1(@testfunction2,a,b,h(3));
[xb_3,yb_3] = bwd1(@testfunction2,a,b,h(3));
[xc1_3,yc1_3] = cen1(@testfunction2,a,b,h(3));
[xc2_3,yc2_3] = cen2(@testfunction2,a,b,h(3));

% Derivatives obtained through each method when h = 0.001

[xf_4,yf_4] = fwd1(@testfunction2,a,b,h(4));
[xb_4,yb_4] = bwd1(@testfunction2,a,b,h(4));
[xc1_4,yc1_4] = cen1(@testfunction2,a,b,h(4));
[xc2_4,yc2_4] = cen2(@testfunction2,a,b,h(4));

% First derivative obtained analytically

x = a:.01:b;
y_prime = 2*x-1;

% Second derivative obtained analytically

y_doubleprime = 2*ones(1,length(x));

% Absolute error for each function for each value of h

% h = 1

errorf_1 = abs(yf_1 - (2*xf_1 - 1));
errorb_1 = abs(yb_1 - (2*xb_1 - 1));
errorc1_1 = abs(yc1_1 -(2*xc1_1 - 1));
errorc2_1 = abs(yc2_1 - 2*ones(1,length(xc2_1)));

% h = 0.1

errorf_2 = abs(yf_2 - (2*xf_2 - 1));
errorb_2 = abs(yb_2 - (2*xb_2 - 1));
errorc1_2 = abs(yc1_2 - (2*xc1_2 - 1));
errorc2_2 = abs(yc2_2 - 2*ones(1,length(xc2_2)));

% h = 0.01

errorf_3 = abs(yf_3 - (2*xf_3 - 1));
errorb_3 = abs(yb_3 - (2*xb_3 - 1));
errorc1_3 = abs(yc1_3 - (2*xc1_3 - 1));
errorc2_3 = abs(yc2_3 - 2*ones(1,length(xc2_3)));

% h = 0.001

errorf_4 = abs(yf_4 - (2*xf_4 - 1));
errorb_4 = abs(yb_4 - (2*xb_4 - 1));
errorc1_4 = abs(yc1_4 - (2*xc1_4 - 1));
errorc2_4 = abs(yc2_4 - 2*ones(1,length(xc2_4)));

% Plot the numerical solution for each step size along with the analytical
% solution

% Forward-Difference

figure(13);
plot(x,y_prime,'r-'); hold on
plot(xf_1,yf_1,'b-')
plot(xf_2,yf_2,'g-')
plot(xf_3,yf_3,'c-')
plot(xf_4,yf_4,'k-')
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [FWD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001','Location','Southeast')
print('P1_13', '-dpng', '-r300');

% Backward-Difference

figure(14);
plot(x,y_prime,'r-'); hold on
plot(xb_1,yb_1,'b-')
plot(xb_2,yb_2,'g-')
plot(xb_3,yb_3,'c-')
plot(xb_4,yb_4,'k-')
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [BWD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001','Location','Southeast')
print('P1_14', '-dpng', '-r300');

% Central-Difference

figure(15);
plot(x,y_prime,'r-'); hold on
plot(xc1_1,yc1_1,'b-')
plot(xc1_2,yc1_2,'g-')
plot(xc1_3,yc1_3,'c-')
plot(xc1_4,yc1_4,'k-')
xlabel('x')
ylabel('f''(x)')
title('Analytical solution vs Numerical solution ( for different values of h) [CD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001','Location','Southeast')
print('P1_15', '-dpng', '-r300');

% Central-Difference (Second derivative)

figure(16);
plot(x,y_doubleprime,'r-'); hold on
plot(xc2_1,yc2_1,'b-')
plot(xc2_2,yc2_2,'g-')
plot(xc2_3,yc2_3,'c-')
plot(xc2_4,yc2_4,'k-')
xlabel('x')
ylabel('f''''(x)')
ylim([1 3])
title('Analytical solution vs Numerical solution ( for different values of h) [CD method]')
legend('Analytical','h = 1','h = 0.1','h = 0.01','h = 0.001','Location','Southeast')
print('P1_16', '-dpng', '-r300');


% Plot of the derivative at x = 2 obtained from each method vs h

% Forward-difference

figure(17);
loglog([h(1),h(2),h(3),h(4)],[errorf_1(3),errorf_2(21),errorf_3(201),...
    errorf_4(2001)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [FWD method]')
print('P1_17', '-dpng', '-r300');

% Backward-difference

figure(18);
loglog([h(1),h(2),h(3),h(4)],[errorb_1(2),errorb_2(20),errorb_3(200),...
    errorb_4(2000)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [BWD method]')
print('P1_18', '-dpng', '-r300');

% Central-difference

figure(19);
semilogx([h(1),h(2),h(3),h(4)],[errorc1_1(2),errorc1_2(20),errorc1_3(200),...
    errorc1_4(2000)])
xlabel('h')
ylabel('absolute error ')
title('Absolute error of f''(x) at x = 2 vs step size (h)  [CD method]')
print('P1_19', '-dpng', '-r300');

% Central-difference (Second Derivative)

figure(20);
semilogx([h(1),h(2),h(3),h(4)],[errorc2_1(2),errorc2_2(20),errorc2_3(200),...
    errorc2_4(2000)])
xlabel('h')
ylabel('absolute error')
title('Absolute error of f''''(x) at x = 2 vs step size (h)  [CD method]')
print('P1_20', '-dpng', '-r300');

% Plot of absolute error on the range

% Forward-difference

figure(21);
semilogy(xf_1,errorf_1,'r-',xf_2,errorf_2,'b-',xf_3,errorf_3,'g-',...
    xf_4,errorf_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [FWD method]')
ylim([10^(-4) 10])
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_21', '-dpng', '-r300');

% Backward-difference

figure(22);
semilogy(xb_1,errorb_1,'r-',xb_2,errorb_2,'b-',xb_3,errorb_3,'g-',...
    xb_4,errorb_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [BWD method]')
ylim([10^(-4) 10])
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_22', '-dpng', '-r300');

% Central-difference

figure(23);
semilogy(xc1_1,errorc1_1,'r-',xc1_2,errorc1_2,'b-',xc1_3,errorc1_3,'g-',...
    xc1_4,errorc1_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''(x) [CD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_23', '-dpng', '-r300');

% Central-difference (Second Derivative)

figure(24);
semilogy(xc2_1,errorc2_1,'r-',xc2_2,errorc2_2,'b-',xc2_3,errorc2_3,'g-',...
    xc2_4,errorc2_4,'k-')
xlabel('x')
ylabel('absolute error')
title('Absolute error of f''''(x) [CD method]')
legend('h = 1','h = 0.1','h = 0.01','h = 0.001')
print('P1_24', '-dpng', '-r300');





