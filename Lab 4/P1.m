% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to apply Crank-Nicolson method to a simple heat
% equation problem

% External m.files used in this code:
% -> ana_solution_P1.m

clear all; clc;

% Mesh parameters

h = 0.01;    % space step size  
k = 0.01;  % time step

x0 = 0;   xfinal = 1;
t0 = 0;   tfinal = 1;

x = x0:h:xfinal ;
t = t0:k:tfinal ;

N = length(x); 
M = length(t); 

% Boundary conditions and Initial condition

w = zeros(N,M);   % w is the approximation of T

w(:,1) = ones(N,1);    % T = 1 , t = 0 and 0 < x < 1
w(1,:) = 0;            % T = 0, t > 0 and x = 0
w(N,:) = 0;            % T = 0, t > 0 and x = 1


% Applying the Crank-Nicolson scheme to the heat equation, the following
% system of linear equations is obtained:
%
%                   Awj+1 = Bwj + C
%

% Building A and B matrices

A_diag = (1+k/h^2)*ones(N-2,1);

A_off_diag = -k/(2*h^2)*ones(N-3,1);

A = diag(A_diag) + diag(A_off_diag,1) + diag(A_off_diag,-1);

B_diag = (1-k/h^2)*ones(N-2,1);

B_off_diag = k/(2*h^2)*ones(N-3,1);

B = diag(B_diag) + diag(B_off_diag,1) + diag(B_off_diag,-1);

C = zeros(N-2,1);

for j = 1:M-2
    
    C(1) = k/h^2*(w(1,j) + w(1,j+1));
    C(N-2) = k/h^2*(w(N,j) + w(N,j+1));
    w(2:N-1,j+1) = A\(B*w(2:N-1,j)+C);
    
end

% Requested plot

figure(1); 
plot(x,w(:,1+1/k),'r*-',x,ana_solution_P1(x,1),'b-','LineWidth',1.5)
grid on;
set(gca,'FontSize',17);
xlabel('x');
ylabel('T(x,1)');
legend('Numerical solution','Analytical Solution')
title('Numerical vs analytical solution at t = 1');
print('P1_1', '-dpng', '-r300');

%%

% Testing of stability of the method

clear all; clc;

% Mesh parameters

h = 0.01;    % space step size  
k = [0.1 0.01 0.001 ];  % time step

colours = ['r-','g-','b-','m-'];
figure(2) ; grid on; hold on;

for q = 1:length(k)

    x0 = 0;   xfinal = 1;
    t0 = 0;   tfinal = 1;

    x = x0:h:xfinal ;
    t = t0:k(q):tfinal ;

    N = length(x); 
    M = length(t); 

    % Boundary conditions and Initial condition

    w = zeros(N,M);   % w is the approximation of T

    w(:,1) = ones(N,1);    % T = 1 , t = 0 and 0 < x < 1
    w(1,:) = 0;            % T = 0, t > 0 and x = 0
    w(N,:) = 0;            % T = 0, t > 0 and x = 1


    % Applying the Crank-Nicolson scheme to the heat equation, the following
    % system of linear equations is obtained:
    %
    %                   Awj+1 = Bwj + C
    %

    % Building A and B matrices

    A_diag = (1+k(q)/h^2)*ones(N-2,1);

    A_off_diag = -k(q)/(2*h^2)*ones(N-3,1);

    A = diag(A_diag) + diag(A_off_diag,1) + diag(A_off_diag,-1);

    B_diag = (1-k(q)/h^2)*ones(N-2,1);

    B_off_diag = k(q)/(2*h^2)*ones(N-3,1);

    B = diag(B_diag) + diag(B_off_diag,1) + diag(B_off_diag,-1);

    C = zeros(N-2,1);

    for j = 1:M-2

        C(1) = k(q)/h^2*(w(1,j) + w(1,j+1));
        C(N-2) = k(q)/h^2*(w(N,j) + w(N,j+1));
        w(2:N-1,j+1) = A\(B*w(2:N-1,j)+C);

    end
    
    % Plots
    
    plot(x,w(:,1+0.1/k(q)),colours(q),'LineWidth',1.5)    
    
end
 
set(gca,'FontSize',17);
xlabel('x')
ylabel('T(x,0.1)')
legend('k = 0.1', 'k = 0.01','k = 0.001','Location','SouthEast')
title('Numerical solution for a range of k''s and fixed h (h=0.01)')
print('P1_2', '-dpng', '-r300');
hold off
    




