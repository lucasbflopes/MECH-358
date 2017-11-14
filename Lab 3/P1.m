% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve the heat equation with precribed boundary
% conditions using finite differences (Forward time centered space scheme)

% External m.files used in this code:
% -> ana_solution.m

clear all;  clc;

% Mesh parameters

h = 0.1;    % space step size  
k = 0.001;  % time step

x0 = 0;   xfinal = 1;
t0 = 0;   tfinal = 1;

x = x0:h:xfinal ;
t = t0:k:tfinal ;

N = length(x); 
M = length(t); 

% Boundary conditions

T0 = 0;
Tfinal = 0;

% Setting Tij

T_x0 = @(x) x.*(1-x);
T = zeros(N,M);
T(:,1) = T_x0(x)';           % Initial condition
T(1,:) = T0*ones(1,M);       % BC @ x = 0   
T(N,:) = Tfinal*ones(1,M);   % BC @ x = 1   

% Computing Tij+1

for j = 1:M-1
    for i = 2:N-1
        
        T(i,j+1) = k/h^2*T(i+1,j) + (1-2*k/h^2)*T(i,j) + k/h^2*T(i-1,j);
        
    end    
end    

% Plot of requested graphs

figure(1);
hold on
plot(x,ana_solution(x,0),'b-',x,ana_solution(x,0.05),'r-',...
    x,ana_solution(x,0.1),'g-',x,ana_solution(x,0.2),'m-',...
    x,ana_solution(x,1),'c-')
plot(x,T(:,1),'b*',x,T(:,1+0.05/k),'r*',x,T(:,1+0.1/k),'g*',x,T(:,1+0.2/k)...
    ,'m*',x,T(:,1+1/k),'c*')
hold off
xlabel('x')
ylabel('T(x)')
legend('t = 0','t = 0.05','t = 0.1','t = 0.2','t = 1')
title('Numerical solution(*) vs Analytical(-) solution of the heat equation')
print('P1_1.png', '-dpng', '-r300');

% Test of instability

% Based on the equation of stability developed in class, for h = 0.1 the
% threshold value of k should be 0.005 
%
%                   k > 0.005   (instability)
%

clear all;  clc

% Mesh parameters

h = 0.1;    % space step size  
k = 0.0051;  % time step

x0 = 0;   xfinal = 1;
t0 = 0;   tfinal = 1;

x = x0:h:xfinal ;
t = t0:k:tfinal ;

N = length(x); 
M = length(t); 

% Boundary conditions

T0 = 0;
Tfinal = 0;

% Setting Tij

T_x0 = @(x) x.*(1-x);
T = zeros(N,M);
T(:,1) = T_x0(x)';           % Initial condition
T(1,:) = T0*ones(1,M);       % BC @ x = 0   
T(N,:) = Tfinal*ones(1,M);   % BC @ x = 1   

% Computing Tij+1

for j = 1:M-1
    for i = 2:N-1
        
        T(i,j+1) = k/h^2*T(i+1,j) + (1-2*k/h^2)*T(i,j) + k/h^2*T(i-1,j);
        
    end    
end    

figure(2);
hold on
plot(x,ana_solution(x,0),'r-',x,ana_solution(x,k*(round(M/10)-1)),...
    x,ana_solution(x,1),'m-')
plot(x,T(:,1),'r*',x,T(:,round(M/10)),'b*',x,T(:,M),'m*')
hold off
xlabel('x')
ylabel('T(x)')
legend('t = 0',['t = ',num2str(k*(round(M/10)-1))],'t = 1')
title(['Numerical solution(*) vs Analytical solution(-) when the stability criterion is not satisfied ('...
    'k = ',num2str(k),')'])
print(['P1_k=',num2str(k),'.png'], '-dpng', '-r300');




