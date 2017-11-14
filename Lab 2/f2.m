function [ f_x ] = f2( x )

% This function returns the value of 
%
% f(x) = { -1   , -2 < x < 0 
%        {  1   , 0 <= x < 2 


f_x = zeros(1,length(x));

for i = 1:length(x)

    if x(i)<0
        f_x(i) = -1;
    else
        f_x(i) = 1;
    end
    
end

end

