% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve Laplace's equation using finite-difference and
% Gauss-Seidel

% External m.files used in this code:
% -> ana_solution_P2.m
% -> RHS_BC.m

%%

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%           Part 1: h = 0.2
%--------------------------------------------------------------------
%--------------------------------------------------------------------

clear all; clc;

% Mesh parameters

h = 0.2;  % space step size both in x and y

x0 = 0;  xfinal = 1;

y0 = 0;  yfinal = 1;

x = x0:h:xfinal;

y = y0:h:yfinal;

N = length(x); % which is also equals to length(y)

% Defining T

T = zeros(N,N);

% Setting first value of Tij ( initial guess) : Tij = 0

T(2:N-1,2:N-1) = 0;

% Setting boundary conditions

T(1,:) = 0;           % T(x=0,y) = 0
T(N,:) = RHS_BC(y);   % T(x=1,y) = f(y)
T(:,1) = 0;           % T(x,y=0) = 0
T(:,N) = 0;           % T(x,y=1) = 0


% Starting method ( from bottom left to top right of the square )

iter = 7;

for k = 1:iter
    for j = 2:N-1
        for i = 2:N-1
            T(i,j) = 1/4*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
        end
    end
end

% Setting analytical solution

analytical_sol = zeros(N,N);

for i = 1:N
    for j = 1:N
        analytical_sol(i,j) = ana_solution_P2(x(i),y(j));
    end
end

error = abs ( T - analytical_sol );

zlevs = 0.1:0.1:.9;  % the contour plot will highlight these level curves

figure(1);
contour(x,y,T',zlevs,'LineWidth',2)
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the numerical solution (h = 0.2)');
print('P2_1', '-dpng', '-r300');

figure(2);
contour(x,y,analytical_sol',zlevs,'LineWidth',2)
set(gca,'FontSize',17);
colormap jet
xlabel('x');
ylabel('y');
title('Contour of the analytical solution (h = 0.2)');
print('P2_2', '-dpng', '-r300');

figure(3);
contour(x,y,error',10,'LineWidth',2)
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the error (h = 0.2)');
print('P2_3', '-dpng', '-r300');

%%

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%           Part 2: h = 0.05
%--------------------------------------------------------------------
%--------------------------------------------------------------------

clear all; clc;

% Mesh parameters

h = 0.05;  % space step size both in x and y

x0 = 0;  xfinal = 1;

y0 = 0;  yfinal = 1;

x = x0:h:xfinal;

y = y0:h:yfinal;

N = length(x); % which is also equals to length(y)

% Defining T

T = zeros(N,N);

% Setting first value of Tij ( initial guess) : Tij = 0

T(2:N-1,2:N-1) = 0;

% Setting boundary conditions

T(1,:) = 0;           % T(x=0,y) = 0
T(N,:) = RHS_BC(y);   % T(x=1,y) = f(y)
T(:,1) = 0;           % T(x,y=0) = 0
T(:,N) = 0;           % T(x,y=1) = 0


% Starting method ( from bottom left to top right of the square )

iter = 7;

for k = 1:iter
    for j = 2:N-1
        for i = 2:N-1
            T(i,j) = 1/4*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
        end
    end
end

% Setting analytical solution

analytical_sol = zeros(N,N);

for i = 1:N
    for j = 1:N
        analytical_sol(i,j) = ana_solution_P2(x(i),y(j));
    end
end

error = abs ( T - analytical_sol );

zlevs = 0.1:0.1:.9;  % the contour plot will highlight these level curves

figure(4);
contour(x,y,T',zlevs,'LineWidth',2)
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the numerical solution (h = 0.05)');
print('P2_4', '-dpng', '-r300');

figure(5);
contour(x,y,analytical_sol',zlevs,'LineWidth',2)
set(gca,'FontSize',17);
colormap jet
xlabel('x');
ylabel('y');
title('Contour of the analytical solution (h = 0.05)');
print('P2_5', '-dpng', '-r300');

figure(6);
contour(x,y,error',10,'LineWidth',2)
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the error (h = 0.05)');
print('P2_6', '-dpng', '-r300');

%%

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%           Part 3: varying the number of iterations
%--------------------------------------------------------------------
%--------------------------------------------------------------------

clear all; clc;

% Mesh parameters

h = 0.05;  % space step size both in x and y

x0 = 0;  xfinal = 1;

y0 = 0;  yfinal = 1;

x = x0:h:xfinal;

y = y0:h:yfinal;

N = length(x); % which is also equals to length(y)

iterations = [10  100 1000];

% Defining T

T = zeros(N,N);

% Setting first value of Tij ( initial guess) : Tij = 0

T(2:N-1,2:N-1) = 0;

% Setting boundary conditions

T(1,:) = 0;           % T(x=0,y) = 0
T(N,:) = RHS_BC(y);   % T(x=1,y) = f(y)
T(:,1) = 0;           % T(x,y=0) = 0
T(:,N) = 0;           % T(x,y=1) = 0


% Starting method ( from bottom left to top right of the square )

cont = 1;

for q = 1:length(iterations)

    for k = 1:iterations(q) 
        for j = 2:N-1
            for i = 2:N-1
                T(i,j) = 1/4*(T(i+1,j) + T(i-1,j) + T(i,j+1) + T(i,j-1));
            end
        end
    end

    % Setting analytical solution

    analytical_sol = zeros(N,N);

    for i = 1:N
        for j = 1:N
            analytical_sol(i,j) = ana_solution_P2(x(i),y(j));
        end
    end

    error = abs ( T - analytical_sol );

    zlevs = 0.1:0.1:.9;  % the contour plot will highlight these level curves

    figure(7+2*(cont-1));
    contour(x,y,T',zlevs,'LineWidth',2)
    set(gca,'FontSize',17);
    xlabel('x');
    ylabel('y');
    title(['Contour of the numerical solution (h = 0.05) (iterations = ',...
        num2str(iterations(q)),')']);
    print(['P2_',num2str(7+2*(cont-1))], '-dpng', '-r300');

    figure(8+2*(cont-1));
    contour(x,y,error',10,'LineWidth',2)
    set(gca,'FontSize',17);
    xlabel('x');
    ylabel('y');
    title(['Contour of the error (h = 0.05) (iterations = ',...
        num2str(iterations(q)),')']);
    print(['P2_',num2str(8+2*(cont-1))], '-dpng', '-r300');

    cont = cont+1;

end
