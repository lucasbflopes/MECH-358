function [ A ] = tridiag( a,b,c,N )

% This function creates a tridiagonal matrix N x N with constants elements
% a on the diagonal, b on the  superdiagonal and c on the subdiagonal

A = zeros(N,N);

for i = 1:N
    for j = 1:N
        if i==j
            A(i,j) = a;
        elseif j == i+1
            A(i,j) = b;
        elseif j == i-1
            A(i,j) = c;
        end
    end
end

end

