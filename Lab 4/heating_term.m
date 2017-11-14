function [ sol ] = heating_term( x,y )

% This function is the heating term function on Problem 3 item c

if 0.3<=x && x<=0.7 && 0.3<=y && y<= 0.7
    
    sol = -0.5;
    
else
    
    sol = 0;
    
end

end

