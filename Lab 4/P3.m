% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to solve Poisson equations with Neumann boundary
% conditions, using finite differences and Gauss-Seidel

% External m.files used in this code:
% -> heating_term.m

%%

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%                    Part 1: Itens A and B
%--------------------------------------------------------------------
%--------------------------------------------------------------------

clear all; clc;

% Mesh parameters

h = 0.1;  % space step size both in x and y

x0 = 0;  xfinal = 1;

y0 = 0;  yfinal = 1;

x = x0:h:xfinal;

y = y0:h:yfinal;

N = length(x); % which is also equals to length(y)

% Defining f(x,y)

f = 1 ;

% Defining phi

phi = zeros(N,N);  % first guess

% Setting Dirichlet boundary condition

phi(N,:) = zeros(1,N);  % phi(1,y) = 0

% Starting method 

for k = 1:700 % 700 iterations
        for j = 1:N
            for i = 1:N-1
                if i == 1 && j == 1

                    phi(i,j) = 1/4*(2*phi(2,1) + 2*phi(1,2)) - ...
                        h^2/4*f; 

                elseif i == 1 && j == N

                    phi(i,j) = 1/4*(2*phi(2,N) + 2*phi(1,N-1)) - ...
                        h^2/4*f; 

                elseif i == 1

                    phi(i,j) = 1/4*(2*phi(2,j) + phi(1,j+1) + phi(1,j-1))...
                         - h^2/4*f; 

                elseif j == 1

                    phi(i,j) = 1/4*(2*phi(i,2) + phi(i+1,1) + phi(i-1,1))...
                         - h^2/4*f;                

                elseif j == N

                    phi(i,j) = 1/4*(2*phi(i,N-1) + phi(i+1,N) + phi(i-1,N))...
                         - h^2/4*f;           

                else

                    phi(i,j) = 1/4*(phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + ...
                        phi(i,j-1)) - h^2/4*f;

                end

            end
        end
end

% Analytical solution

g = @(x,y) 1/2*(x.^2-1);

analytical_sol = zeros(N,N);

for i=1:N
    for j=1:N
        analytical_sol(i,j) = g(x(i),y(j));
    end
end

% Contour plots

zlevels = -[ .05 .1 .15 .2 .25 .3 .35 .4 .45];

figure(1);
contour(x,y,phi',zlevels,'LineWidth',2,'ShowText','on')
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the numerical solution (h = 0.1)');
print('P3_1', '-dpng', '-r300');

figure(2);
contour(x,y,analytical_sol',zlevels,'LineWidth',2,'ShowText','on')
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title('Contour of the analytical solution');
print('P3_2', '-dpng', '-r300');



%%

%--------------------------------------------------------------------
% Computing the error @ (x,y)=(0.5,0.6) for different # of iterations
%--------------------------------------------------------------------

clear all;  clc;

% Mesh parameters

h = 0.1;  % space step size both in x and y

x0 = 0;  xfinal = 1;

y0 = 0;  yfinal = 1;

x = x0:h:xfinal;

y = y0:h:yfinal;

N = length(x); % which is also equals to length(y)

% Defining phi

phi = zeros(N,N);

% Defining f(x,y)

f = @(x,y) 1 ;

% Setting Dirichlet boundary condition

phi(N,:) = zeros(1,N);  % phi(1,y) = 0

% Starting method 

iter = [10 50 100 500 1000 5000 10000];

error = zeros(1,length(iter));

for q = 1:length(iter)
    for k = 1:iter(q) % iterations
        for j = 1:N
            for i = 1:N-1
                if i == 1 && j == 1

                    phi(i,j) = 1/4*(2*phi(2,1) + 2*phi(1,2)) - ...
                        h^2/4*f(x(1),y(1)); 

                elseif i == 1 && j == N

                    phi(i,j) = 1/4*(2*phi(2,N) + 2*phi(1,N-1)) - ...
                        h^2/4*f(x(1),y(N)); 

                elseif i == 1

                    phi(i,j) = 1/4*(2*phi(2,j) + phi(1,j+1) + phi(1,j-1))...
                         - h^2/4*f(x(1),y(j)); 

                elseif j == 1

                    phi(i,j) = 1/4*(2*phi(i,2) + phi(i+1,1) + phi(i-1,1))...
                         - h^2/4*f(x(i),y(1));                

                elseif j == N

                    phi(i,j) = 1/4*(2*phi(i,N-1) + phi(i+1,N) + phi(i-1,N))...
                         - h^2/4*f(x(i),y(N));           

                else

                    phi(i,j) = 1/4*(phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + ...
                        phi(i,j-1)) - h^2/4*f(x(i),y(j));

                end

            end
        end
    end

    % Analytical solution

    g = @(x,y) 1/2*(x.^2-1);

    analytical_sol = zeros(N,N);

    for i=1:N
        for j=1:N
            analytical_sol(i,j) = g(x(i),y(j));
        end
    end

    error(q) = abs(analytical_sol(6,7)-phi(6,7));

end

figure(3);
loglog(iter,error,'LineWidth',1.5); grid on;
set(gca,'FontSize',17);
xlabel('# of iterations')
ylabel('absolute error')
title('Absolute error of numerical solution at (x,y) = (0.5,0.6)')
print('P3_3', '-dpng', '-r300');

%% 

%--------------------------------------------------------------------
%--------------------------------------------------------------------
%                    Part 2: Itens C and D
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

% Defining phi

phi = zeros(N,N);

% Defining f(x,y)

% f(x,y) is the auxiliary function heat_term.m

% Setting Dirichlet boundary condition

phi(N,:) = zeros(1,N);  % phi(1,y) = 0

% Starting method 

for k = 1:12000 % 12000 iterations
    for j = 1:N
        for i = 1:N-1
            if i == 1 && j == 1
                
                phi(i,j) = 1/4*(2*phi(2,1) + 2*phi(1,2)) - ...
                    h^2/4*heating_term(x(1),y(1)); 
                
            elseif i == 1 && j == N
                
                phi(i,j) = 1/4*(2*phi(2,N) + 2*phi(1,N-1)) - ...
                    h^2/4*heating_term(x(1),y(N)); 
                
            elseif i == 1
                
                phi(i,j) = 1/4*(2*phi(2,j) + phi(1,j+1) + phi(1,j-1))...
                     - h^2/4*heating_term(x(1),y(j)); 
                 
            elseif j == 1
                
                phi(i,j) = 1/4*(2*phi(i,2) + phi(i+1,1) + phi(i-1,1))...
                     - h^2/4*heating_term(x(i),y(1));                
                 
            elseif j == N
                
                phi(i,j) = 1/4*(2*phi(i,N-1) + phi(i+1,N) + phi(i-1,N))...
                     - h^2/4*heating_term(x(i),y(N));           
           
            else
                
                phi(i,j) = 1/4*(phi(i+1,j) + phi(i-1,j) + phi(i,j+1) + ...
                    phi(i,j-1)) - h^2/4*heating_term(x(i),y(j));
                
            end
            
        end
    end
end

% Contour plot

zlevels = [ .01 .02 .03 .04 .045 .048];

figure(4);
contour(x,y,phi',zlevels,'LineWidth',2,'ShowText','on')
set(gca,'FontSize',17);
xlabel('x');
ylabel('y');
title(['Contour of the numerical solution (h = ',num2str(h),')']);
print('P3_4', '-dpng', '-r300');




