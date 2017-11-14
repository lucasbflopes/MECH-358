function [ sol ] = RHS_BC( T )

% This function applies the following function to each component of the
% array T
%                         { 0 , 0 <= y <= 1/2
%                 f(y) =  {
%                         { 1 , 1/2< y < 1          
%

sol = zeros(length(T),1);

for i = 1:length(T)

    if T(i)>=0 && T(i)<=1/2
        sol(i) = 0;
    elseif T(i)>1/2 && T(i)<=1
        sol(i) = 1;
    else
        sol(i) = 0;
    end

end

end

