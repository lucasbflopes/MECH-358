% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to study the convergence of Fourier series for a
% triangular and square wave functions

% External m.files used in this code:
% -> S1.m
% -> S2.m
% -> f2.m

clear all ; clc

x = -2:.1:2;

% Function f1 : f(x) = |x|   , -2<=x<=2

f_1 = abs(x);

% Function f2: f(x) = { -1 on  -2<=x<=0        
%                     {  1 on   0<=x<=2  

f_2 = f2(x);

% Given a value of N, we calculate the partial sum of the Fourier Serie of
% each function up to the Nth term

% N = 5

f1_5 = S1(x,5);
f2_5 = S2(x,5);

% N = 10

f1_10 = S1(x,10);
f2_10 = S2(x,10);

% N = 30

f1_30 = S1(x,30);
f2_30 = S2(x,30);

% N = 100

f1_100 = S1(x,100);
f2_100 = S2(x,100);

% Plot the Fourier Series of f1(x) for different values of N

figure(1);
plot(x,f1_5,'r',x,f1_10,'b',x,f1_30,'g',x,f1_100,'k')
xlabel('x')
ylabel('f(x)')
legend('N = 5','N = 10','N = 30','N = 100')
title('Fourier Series of f1(x) for different values of N')
print('P3_1', '-dpng', '-r300');

% Plot the Fourier Series of f2(x) for different values of N

hold on
figure(2);
plot(x,f2_5,'r',x,f2_10,'b',x,f2_30,'g',x,f2_100,'k')
xlabel('x')
ylabel('f(x)')
title('Fourier Series of f2(x) for different values of N')
legend('N = 5','N = 10','N = 30','N = 100','Location','Northwest')
print('P3_2', '-dpng', '-r300');

% Computing the relative error of the series for each N

% N = 5

error_1_5 = abs(S1(x,5)-abs(x));
error_2_5 = abs(S2(x,5)-f2(x));

% N = 10

error_1_10 = abs(S1(x,10)-abs(x));
error_2_10 = abs(S2(x,10)-f2(x));

% N = 30

error_1_30 = abs(S1(x,30)-abs(x));
error_2_30 = abs(S2(x,30)-f2(x));

% N = 100

error_1_100 = abs(S1(x,100)-abs(x));
error_2_100 = abs(S2(x,100)-f2(x));

% Plot of the absolute error for both functions

figure(3)
plot(x,error_1_5,'r-',x,error_1_10,'b-',x,error_1_30,'g-',x,error_1_100,'k')
xlabel('x')
ylabel('absolute error')
legend('N = 5','N = 10','N = 30','N = 100')
title('Absolute error of f1(x) for different values of N')
print('P3_3', '-dpng', '-r300');


figure(4)
plot(x,error_2_5,'r-',x,error_2_10,'b-',x,error_2_30,'g-',x,error_2_100,'k')
xlabel('x')
ylabel('absolute error')
legend('N = 5','N = 10','N = 30','N = 100')
title('Absolute error of f2(x) for different values of N')
print('P3_4', '-dpng', '-r300');

% Computing the relative error of S_n(x) @ x = 1/2 for the different 
% values of N (results in %)

% OBS: f1(x = 1/2) = abs(1/2) = 1/2
%      f2(x = 1/2) = 1 

% N = 5

rel_error_1_5 = (S1(1/2,5)-1/2)/(1/2)*100;
rel_error_2_5 = (S2(1/2,5)-1)*100;

% N = 10

rel_error_1_10 = (S1(1/2,10)-1/2)/(1/2)*100;
rel_error_2_10 = (S2(1/2,10)-1)*100;


% N = 30

rel_error_1_30 = (S1(1/2,30)-1/2)/(1/2)*100;
rel_error_2_30 = (S2(1/2,30)-1)*100;


% N = 100

rel_error_1_100 = (S1(1/2,100)-1/2)/(1/2)*100;
rel_error_2_100 = (S2(1/2,100)-1)*100;

% Table of values

disp('f1(x)')

disp(['N = 5   || Sn(x = 1/2) =',num2str(S1(1/2,5)),...
    ' || relative error =',num2str(rel_error_1_5)])

disp(['N = 10  || Sn(x = 1/2) =',num2str(S1(1/2,10)),...
    ' || relative error =',num2str(rel_error_1_10)])

disp(['N = 30  || Sn(x = 1/2) =',num2str(S1(1/2,30)),...
    ' || relative error =',num2str(rel_error_1_30)])

disp(['N = 100 || Sn(x = 1/2) =',num2str(S1(1/2,100)),...
    ' || relative error =',num2str(rel_error_1_100)])

disp('f2(x)')

disp(['N = 5   || Sn(x = 1/2) =',num2str(S2(1/2,5)),...
    ' || relative error =',num2str(rel_error_2_5)])

disp(['N = 10  || Sn(x = 1/2) =',num2str(S2(1/2,10)),...
    ' || relative error =',num2str(rel_error_2_10)])

disp(['N = 30  || Sn(x = 1/2) =',num2str(S2(1/2,30)),...
    ' || relative error =',num2str(rel_error_2_30)])

disp(['N = 100 || Sn(x = 1/2) =',num2str(S2(1/2,100)),...
    ' || relative error =',num2str(rel_error_2_100)])






