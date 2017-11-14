% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to provide a deepen knowledge about the study of
% truncation error and the dependence of the truncation error on the mesh
% size

clear all;  clc;

% --------------------------------------------------------------
% ------------------- Item a  [ O(h2) scheme]-------------------
% --------------------------------------------------------------

h = [ 0.1 0.05 0.01 0.005];   % space step size

error = zeros(length(h),2); % Used later to store the relative error @ y(0.5) 
%                             for different h's and both schemes

for i = 1:length(h)
    
    % Mesh parameters

    x0 = 0 ; xfinal = 1 ;

    x = x0:h(i):xfinal;

    N = length(x);

    % Boundary conditions

    y = zeros(N,1);
    y(1) = 0;     % BC y(0) = 0
    y(N) = 0;     % BC y(1) = 0      


    % By using the central difference scheme, the ode becomes:
    %
    %           y(xi-1) - 2*y(xi) + y(xi+1) = -h^2*exp(xi)
    %
    %                   A*y = b

    % Building A matrix and b vector

    diagonal = -2*ones(N-2,1);
    upperdiagonal = 1*ones(N-3,1);
    lowerdiagonal = 1*ones(N-3,1);

    A = diag(diagonal) + diag(upperdiagonal,1) + diag(lowerdiagonal,-1);

    b = zeros(N-2,1);
    b(1) = -h(i)^2*exp(x(1)) - y(1);
    b(N-2) = -h(i)^2*exp(x(N)) - y(N);
    b(2:N-3) = -h(i)^2*exp(x(2:N-3));

    % Solving the system

    y(2:N-1) = A\b ;

    % Storing the value at the midpoint

    error(i,1) = abs(0.2104196 - y(0.5/h(i) + 1))/0.2104196;

end

% --------------------------------------------------------------
% ------------------- Item b  [ O(h4) scheme]-------------------
% --------------------------------------------------------------

for i = 1:length(h)

    % Mesh parameters

    x0 = 0 ; xfinal = 1 ;

    x = x0:h(i):xfinal;

    N = length(x);

    % Boundary conditions

    y = zeros(N,1);
    y(1) = 0;     % BC y(0) = 0
    y(N) = 0;     % BC y(1) = 0   

    % By using the the 4th order scheme on HW2 P2, the ode becomes:
    %
    %  y(xi-1)-2*y(xi)+y(xi+1) = -h^2/12*(f(xi+1)+10*f(xi)+f(xi-1)
    %
    %                   A*y = B*f + c

    % Building A and B matrices as well as f and c vectors

    diagonal_a = -2*ones(N-2,1);
    upperdiagonal_a = 1*ones(N-3,1);
    lowerdiagonal_a = 1*ones(N-3,1);

    diagonal_b = -10*h(i)^2/12*ones(N-2,1);
    upperdiagonal_b = -h(i)^2/12*ones(N-3,1);
    lowerdiagonal_b = -h(i)^2/12*ones(N-3,1);

    A = diag(diagonal_a) + diag(upperdiagonal_a,1) + diag(lowerdiagonal_a,-1);
    B = diag(diagonal_b) + diag(upperdiagonal_b,1) + diag(lowerdiagonal_b,-1);

    f = exp(x)';

    c = zeros(N-2,1);
    c(1) =  - y(1) -h(i)^2/12*f(1);
    c(N-2) = - y(N) -h(i)^2/12*f(N);

    % Solving the system

    y(2:N-1) = A\(B*f(2:N-1)+c) ;

    % Storing the relative error at the midpoint

    error(i,2) = abs(0.2104196 - y(0.5/h(i) + 1))/0.2104196;

end

% Plotting relative error @ x = 0.5 vs step size for each scheme

figure(1);
loglog(h,error(:,1),'r*-',h,error(:,2),'bo-')
xlabel('h')
ylabel('relative error')
axis([0.5*10^(-2) 10^(-1) 10^(-8) 10^(-1)])
legend('O(h2) scheme','O(h4) scheme','Location','Northwest')
title('Relative error at x = 0.5 as a function of step size(h)')
print('P3_1', '-dpng', '-r300');

% Table of associated error and step size

display(['h',' relative error (O(h^2)scheme)',' relative error (O(h^4)scheme)'])
display([h',error(:,1),error(:,2)]) 

% ----------------------------------------------
% ------------------- Item c -------------------
% ----------------------------------------------

clear all; 

h = [0.25 0.1 0.05];   % space step size

error_relative = zeros(length(h),1);

figure(2)
hold on                        
graph = {'r*-' 'b+-' 'k.-'};

for i = 1:length(h)

    % Mesh parameters

    x0 = 0 ; xfinal = 1 ;

    x = x0:h(i):xfinal;

    N = length(x);

    % Boundary conditions

    y = zeros(N,1);
    y(1) = 0;     % BC y(0) = 0
    y(N) = 0;     % BC y(1) = 0   

    % By using the the 4th order scheme on HW2 P2, the ode becomes:
    %
    %  y(xi-1)-2*y(xi)+y(xi+1) = -h^2/12*(f(xi+1)+10*f(xi)+f(xi-1)
    %
    %                   A*y = B*f + c

    % Building A and B matrices as well as f and c vectors

    diagonal_a = -2*ones(N-2,1);
    upperdiagonal_a = 1*ones(N-3,1);
    lowerdiagonal_a = 1*ones(N-3,1);

    diagonal_b = -10*h(i)^2/12*ones(N-2,1);
    upperdiagonal_b = -h(i)^2/12*ones(N-3,1);
    lowerdiagonal_b = -h(i)^2/12*ones(N-3,1);

    A = diag(diagonal_a) + diag(upperdiagonal_a,1) + diag(lowerdiagonal_a,-1);
    B = diag(diagonal_b) + diag(upperdiagonal_b,1) + diag(lowerdiagonal_b,-1);

    f = (sqrt(x).*sin(10*x))';

    c = zeros(N-2,1);
    c(1) =  - y(1) -h(i)^2/12*f(1);
    c(N-2) = - y(N) -h(i)^2/12*f(N);

    % Solving the system

    y(2:N-1) = A\(B*f(2:N-1)+c) ;

    % Storing the relative and absolute errors at the midpoint

    error_relative(i) = abs(-0.00490196320871196 - y(0.5/h(i) + 1))...
        /0.00490196320871196;
    
    % Plotting y(x) for each step size

    plot(x,y,graph{i})
    xlabel('x')
    ylabel('y')

end

hold off
legend('h = 0.25','h = 0.1','h = 0.05','Location','Northwest');
title('Numerical solution for each step size (h)')
print('P3_2', '-dpng', '-r300');


% Plotting the relative error @ x = 0.5 vs step size

figure(3)
semilogy(h,error_relative,'b*-')
xlabel('h')
ylabel('relative error')
title('Relative error @ x = 0.5 as a function of the step size (h)')
print('P3_3', '-dpng', '-r300');

% Table of associated error and step size

display(['h',' relative error'])
display([h',error_relative]) 


