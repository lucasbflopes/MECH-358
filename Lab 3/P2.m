% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to apply finite difference methods for the heat equation
% for more complex boundary conditions

% ------------------------------------------------------
% ----------------------- CASE 1 -----------------------
% ------------------------------------------------------

clear all; clc;

% Mesh parameters

L = 1;

h = 0.1;    % space step size  
k = 0.001;  % time step

if 2 > h^2/k
    display('This method might not be stable. Consider changing your meshing')
    return
end

x0 = 0;   xfinal = L;
t0 = 0;   tfinal = 10;

x = x0:h:xfinal ;
t = t0:k:tfinal ;

N = length(x); 
M = length(t);

% Boundary conditions, Initial condition, Heat source/sink function

F = 0;              % Initial condition  T(x,0)
G = 0;              % Boundary condition T(0,t)
H = 0;              % Boundary condition T(L,t)

S = @(t) sin(t);    % Heat source/sink equation

% Setting Tij  

T = zeros(N,M);
T(:,1) = F*ones(N,1); 
T(1,:) = G*ones(1,M);
T(N,:) = H*ones(1,M);

% Calculating Tij+1

for j = 1:M-1
    for i = 2:N-1
        
        T(i,j+1) = k/h^2*T(i+1,j) + (1-2*k/h^2)*T(i,j) + k/h^2*T(i-1,j)+...
        k*S(j*k);
        
    end    
end     


% Plot of requested graphs

figure(1);
plot(x,T(:,1),'b*-')
xlabel('x')
ylabel('T(x,0)')
title('Numerical solution of the heat equation (t=0)')
print('P2_1', '-dpng', '-r300');

figure(2);
plot(x,T(:,1+0.01/k),'b*-')
xlabel('x')
ylabel('T(x,0.01)')
title('Numerical solution of the heat equation (t=0.01)')
print('P2_2', '-dpng', '-r300');

figure(3);
plot(x,T(:,1+0.03/k),'b*-')
xlabel('x')
ylabel('T(x,0.03)')
title('Numerical solution of the heat equation (t=0.03)')
print('P2_3', '-dpng', '-r300');

figure(4);
plot(x,T(:,1+0.06/k),'b*-')
xlabel('x')
ylabel('T(x,0.06)')
title('Numerical solution of the heat equation (t=0.06)')
print('P2_4', '-dpng', '-r300');

figure(5);
plot(x,T(:,1+1/k),'b*-')
xlabel('x')
ylabel('T(x,1)')
title('Numerical solution of the heat equation (t=1)')
print('P2_5', '-dpng', '-r300');

figure(6);
plot(x,T(:,1+10/k),'b*-')
xlabel('x')
ylabel('T(x,10)')
title('Numerical solution of the heat equation (t=10)') 
print('P2_6', '-dpng', '-r300');

% Analysis of the oscillation of the midpoint

figure(7);
plot(t,T(1+0.5/h,:))
xlabel('t')
ylabel('T(0.5,t)') 
title('Temperature vs time @ x = 0.5')
axis([0 10 -.15 .15])
print('P2_7', '-dpng', '-r300');

max_amplitude = max(T(1+.5/h,:));
min_amplitude = min(T(1+.5/h,:));

display(['The midpoint oscillates between ',num2str(min_amplitude),...
    ' and ',num2str(max_amplitude)])



% ------------------------------------------------------
% ----------------------- CASE 2 -----------------------
% ------------------------------------------------------


clear all ; 

% Mesh parameters

L = 1;

h = 0.1;    % space step size  
k = 0.001;  % time step

if 2 > h^2/k
    display('This method might not be stable. Consider changing your mesh')
    return
end

x0 = 0;   xfinal = L;
t0 = 0;   tfinal = 10;

x = x0:h:xfinal ;
t = t0:k:tfinal ;

N = length(x); 
M = length(t); 

% Boundary conditions, Initial condition, Heat source/sink function

F = @(x) x/L;               % Initial condition  T(x,0)
G = @(t) 0.5*sin(2*pi.*t);  % Boundary condition T(0,t)
H = 1;                      % Boundary condition T(L,t)

S = 0;                      % Heat source/sink equation

% Setting Tij 

T = zeros(N,M);
T(:,1) = F(x)'; 
T(1,:) = G(t);
T(N,:) = H*ones(1,M);

% Calculating Tij+1

for j = 1:M-1
    for i = 2:N-1
        
        T(i,j+1) = k/h^2*T(i+1,j) + (1-2*k/h^2)*T(i,j) + k/h^2*T(i-1,j);
        
    end    
end    

figure(8);
hold on
plot(x,T(:,1),'r*-')
plot(x,T(:,1+0.25/k),'g*-')
plot(x,T(:,1+0.5/k),'b*-')
plot(x,T(:,1+0.75/k),'c*-')
plot(x,T(:,1+1/k),'k*-')
plot(x,T(:,1+10/k),'m*-')
hold off
xlabel('x')
ylabel('T(x)')
legend('t = 0','t = 0.25','t = 0.5','t = 0.75','t = 1','t = 10','Location'...
    ,'Southeast')
title('Numerical solution of the heat equation')
print('P2_8', '-dpng', '-r300');

% Temperature at x = 0.5

figure(9);
plot(t,T(1+0.5/h,:))
xlabel('t')
ylabel('T(0.5,t)') 
title('Temperature vs time @ x = 0.5')
print('P2_9', '-dpng', '-r300');

max_amplitude = max(T(1+.5/h,:));
min_amplitude = min(T(1+.5/h,:));

display(['The midpoint oscillates between ',num2str(min_amplitude),...
    ' and ',num2str(max_amplitude)])


% Looking for the wall thickness which causes dT/dx @ inner wall < 1

clear all; 

for q = 1:1000  % 1000 is an abitrary value
    
    % Mesh parameters

    L = 1+q*0.01;   % sweeping L in order to find the correct one

    h = 0.1;    % space step size  
    k = 0.001;  % time step

    if 2 > h^2/k
        display('This method might not be stable. Consider changing your mesh')
        return
    end

    x0 = 0;   xfinal = L;
    t0 = 0;   tfinal = 10;

    x = x0:h:xfinal ;
    t = t0:k:tfinal ;

    N = length(x); 
    M = length(t); 

    % Boundary conditions, Initial condition, Heat source/sink function

    F = @(x) x/L;               % Initial condition  T(x,0)
    G = @(t) 0.5*sin(2*pi.*t);  % Boundary condition T(0,t)
    H = 1;                      % Boundary condition T(L,t)

    S = 0;                      % Heat source/sink equation

    % Setting Tij 

    T = zeros(N,M);
    T(:,1) = F(x)'; 
    T(1,:) = G(t);
    T(N,:) = H*ones(1,M);

    % Calculating Tij+1

    for j = 1:M-1
        for i = 2:N-1

            T(i,j+1) = k/h^2*T(i+1,j) + (1-2*k/h^2)*T(i,j) + k/h^2*T(i-1,j);

        end    
    end    

    % Derivative (dT/dx) @ x = L (backward-difference)

    D = (T(N,:) - T(N-1,:))/h;

    if max(D)<1
        display(['The wall thickness should be between ',num2str(L-0.01),...
            ' and ',num2str(L)])
        break
    end
end










