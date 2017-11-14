% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to build an improved method that does a better job than
% the Crank-Niconsol method for heat equation problems discontinuities.

% External m.files used in this code:
% -> ana_solution_P1.m

clear all;  clc;

% ----------------------------------------------------------------
% --------------------  ITEM A -----------------------------------
% ----------------------------------------------------------------

% Mesh parameters

x0 = 0 ;      xN = 1 ;
t0 = 0 ;      tM = 1 ;

h = 0.01; %  space step-size
k = 0.01; %  time step-size

x = x0:h:xN ; 
t = t0:k:tM ;

N = length(x) ;
M = length(t) ;

% Setting numerical solution and Initial/Boundary conditions

T = zeros(N,M);

T(:,1) = 1;    % T(x,0) = 1     IC
T(1,:) = 0;    % T(0,t) = 0     BC
T(N,:) = 0;    % T(1,t) = 0     BC


% Building A matrix from the scheme

diag_main = 1+2*k/h^2*ones(1,N-2) ;
diag_off = -k/h^2*ones(1,N-3) ;

A = diag(diag_main,0) + diag(diag_off,1) + diag(diag_off,-1);

w_j_plus = zeros(N-2,1);
w_j = zeros(N-2,1);

for j = 1:M-1
    w_j(2:N-3) = T(3:N-2,j);
    w_j(1) = T(2,j)+k/h^2*T(1,j+1);
    w_j(N-2) = T(N-1,j)+k/h^2*T(N,j+1);
    w_j_plus = A\w_j;
    T(2:N-1,j+1) = w_j_plus;
end

figure(1);
plot(x,T(:,M),'r-',x,ana_solution_P1(x,1),'b-','LineWidth',1.2)
set(gca,'FontSize',12);
grid on;
xlabel('x')
ylabel('T(x,1)')
legend('Numerical','Analytical')
title('Implicit scheme applied to the heat equation (h=k=0.01)')
print('P1_1', '-dpng', '-r300');

%%

% Verifying the order of the method as well as the stability/instability

clear all; clc;

% Mesh parameters

x0 = 0 ;      xN = 1 ;
t0 = 0 ;      tM = 1 ;

h = 0.01; %  space step-size
k = [0.1 0.05 0.01 0.001]; %  time step-size

error = zeros(1,length(k));

figure(2);  r = {'b-','r-','g-','m-'};
hold on;

for i = 1:length(k)

    x = x0:h:xN ; 
    t = t0:k(i):tM ;

    N = length(x) ;
    M = length(t) ;

    % Setting numerical solution and Initial/Boundary conditions

    T = zeros(N,M);

    T(:,1) = 1;    % T(x,0) = 1     IC
    T(1,:) = 0;    % T(0,t) = 0     BC
    T(N,:) = 0;    % T(1,t) = 0     BC

    % Building A matrix from the scheme

    diag_main = 1+2*k(i)/h^2*ones(1,N-2) ;
    diag_off = -k(i)/h^2*ones(1,N-3) ;

    A = diag(diag_main,0) + diag(diag_off,1) + diag(diag_off,-1);

    w_j_plus = zeros(N-2,1);
    w_j = zeros(N-2,1);

    for j = 1:M-1
        w_j(2:N-3) = T(3:N-2,j);
        w_j(1) = T(2,j)+k(i)/h^2*T(1,j+1);
        w_j(N-2) = T(N-1,j)+k(i)/h^2*T(N,j+1);
        w_j_plus = A\w_j;
        T(2:N-1,j+1) = w_j_plus;
    end

    error(i) = abs(T(1+.5/h,M) - ana_solution_P1(.5,1)); % error of T(0.5,1)
   
    plot(x,T(:,M),r{i},'LineWidth',1.2)

end

set(gca,'FontSize',12);
grid on
xlabel('x')
ylabel('T(x,1)')
legend(['k = ',num2str(k(1))],['k = ',num2str(k(2))],...
    ['k = ',num2str(k(3))],['k = ',num2str(k(4))])
title('Numerical solution at t = 1 for different time-steps (h = 0.01)')
print('P1_2', '-dpng', '-r300');
hold off

figure(3);
plot(k,error,'LineWidth',1.2)
set(gca,'FontSize',12);
grid on;
xlabel('k')
ylabel('absolute error')
title('Error at T(0.5,1) as a function of time-step ( h = 0.01 )')
print('P1_3', '-dpng', '-r300');

%%

% ----------------------------------------------------------------
% --------------------  ITEM B -----------------------------------
% ----------------------------------------------------------------

clear all ; clc; 

% Mesh parameters

x0 = 0 ;      xN = 1 ;
t0 = 0 ;      tM = 1 ;

h = 0.05; %  space step-size
k = 0.1; %  time step-size

x = x0:h:xN ; 
t = t0:2*k:tM ;

N = length(x) ;
M = length(t) ;

% Setting numerical solution and Initial/Boundary conditions

T = zeros(N,M);

T(:,1) = 1;    % T(x,0) = 1     IC
T(1,:) = 0;    % T(0,t) = 0     BC
T(N,:) = 0;    % T(1,t) = 0     BC


% Building M1 and M2 matrices from the scheme

diag_main_A = 1+2*k/h^2*ones(1,N-2) ;
diag_off_A = -k/h^2*ones(1,N-3) ;

A = diag(diag_main_A,0) + diag(diag_off_A,1) + diag(diag_off_A,-1);

diag_main_M1 = 1+2*2*k/h^2*ones(1,N-2) ;
diag_off_M1 = -2*k/h^2*ones(1,N-3) ;

M1 = diag(diag_main_M1,0) + diag(diag_off_M1,1) + diag(diag_off_M1,-1);

M2 = A^2;

wj_2plus_2steps = zeros(N-2,1);
wj_2plus_1step = zeros(N-2,1);

wj_1step = zeros(N-2,1);
wj_2step = zeros(N-2,1);

for j = 1:M-1
    
    wj_1step(2:N-3) = T(3:N-2,j);
    wj_1step(1) = T(2,j)+2*k/h^2*T(1,j+1);
    wj_1step(N-2) = T(N-1,j)+2*k/h^2*T(N,j+1);

    wj_2step(2:N-3) = T(3:N-2,j);
    wj_2step(1) = T(2,j)+k/h^2*T(1,j+1);
    wj_2step(N-2) = T(N-1,j)+k/h^2*T(N,j+1);

    wj_2plus_1step = M1\wj_1step;
    wj_2plus_2steps = M2\wj_2step;
    
    T(2:N-1,j+1) = 2*wj_2plus_2steps - wj_2plus_1step;
    
end

figure(4);
plot(x,T(:,M),'r-',x,ana_solution_P1(x,1),'b-','LineWidth',1.2)
set(gca,'FontSize',12);
grid on
xlabel('x')
ylabel('T(x,1)')
legend('Numerical','Analytical')
title('Extrapolation scheme applied to the heat equation (k=0.1,h=0.05)')
print('P1_4', '-dpng', '-r300');

%%

% Verifying the order of the method as well as the stability/instability

clear all; clc;

% Mesh parameters

x0 = 0 ;      xN = 1 ;
t0 = 0 ;      tM = 1 ;

h = 0.01; %  space step-size
k = [0.1 0.05 0.01 0.001]; %  time step-size

error = zeros(1,length(k));

figure(5);  r = {'b-','r-','g-','m-'};
hold on;

for i = 1:length(k)

    x = x0:h:xN ; 
    t = t0:2*k(i):tM ;

    N = length(x) ;
    M = length(t) ;

    % Setting numerical solution and Initial/Boundary conditions

    T = zeros(N,M);

    T(:,1) = 1;    % T(x,0) = 1     IC
    T(1,:) = 0;    % T(0,t) = 0     BC
    T(N,:) = 0;    % T(1,t) = 0     BC

    % Building M1 and M2 matrices from the scheme

    diag_main_A = 1+2*k(i)/h^2*ones(1,N-2) ;
    diag_off_A = -k(i)/h^2*ones(1,N-3) ;

    A = diag(diag_main_A,0) + diag(diag_off_A,1) + diag(diag_off_A,-1);

    diag_main_M1 = 1+2*2*k(i)/h^2*ones(1,N-2) ;
    diag_off_M1 = -2*k(i)/h^2*ones(1,N-3) ;

    M1 = diag(diag_main_M1,0) + diag(diag_off_M1,1) + diag(diag_off_M1,-1);

    M2 = A^2;

    wj_2plus_2steps = zeros(N-2,1);
    wj_2plus_1step = zeros(N-2,1);
    
    wj_1step = zeros(N-2,1);
    wj_2step = zeros(N-2,1);

    for j = 1:M-1

        wj_1step(2:N-3) = T(3:N-2,j);
        wj_1step(1) = T(2,j)+2*k(i)/h^2*T(1,j+1);
        wj_1step(N-2) = T(N-1,j)+2*k(i)/h^2*T(N,j+1);

        wj_2step(2:N-3) = T(3:N-2,j);
        wj_2step(1) = T(2,j)+k(i)/h^2*T(1,j+1);
        wj_2step(N-2) = T(N-1,j)+k(i)/h^2*T(N,j+1);

        wj_2plus_1step = M1\wj_1step;
        wj_2plus_2steps = M2\wj_2step;

        T(2:N-1,j+1) = 2*wj_2plus_2steps - wj_2plus_1step;

    end

    error(i) = abs(T(1+.5/h,M) - ana_solution_P1(.5,1)); % error of T(0.5,1)
   
    plot(x,T(:,M),r{i},'LineWidth',1.2)

end

set(gca,'FontSize',12);
grid on
xlabel('x')
ylabel('T(x,1)')
legend(['k = ',num2str(k(1))],['k = ',num2str(k(2))],...
    ['k = ',num2str(k(3))],['k = ',num2str(k(4))])
title('Numerical solution at t = 1 for different time-steps (h = 0.01)')
print('P1_5', '-dpng', '-r300');
hold off

figure(6);
plot(k,error,'LineWidth',1.2)
set(gca,'FontSize',12);
grid on;
xlabel('k')
ylabel('absolute error')
title('Error at T(0.5,1) as a function of time-step ( h = 0.01 )')
print('P1_6', '-dpng', '-r300');





