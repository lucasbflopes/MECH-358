% Group: Lucas Bragança Fernandes Lopes and Calvin Sy       
% Student numbers: 56412142 and 57604100

% This code aims to decode a secret message

% External m.files used in this code:
% -> decipher.m
% -> zerobess.m

% URL of the webpage on where the code that finds the zeroes of the bessel
% function were downloaded:

% http://www.mathworks.com/matlabcentral/fileexchange/26639-zerobess

clear all; clc;

% Importing data

sampled_values = csvread('secret.csv.txt');

x = sampled_values(:,1) ;
y = sampled_values(:,2) ;

M = length(x);
h = 0.0005;  % x is sampled at 0.0005 intervals 

% First, assume we have 500 terms. Therefore we need the first 500
% zeroes of the first-kind zero-order bessel function.

zeroes = zerobess('J',0,500);

N = length(zeroes);

% Simpson's rule

Ck = zeros(1,N);

for i = 1:length(zeroes)
    
    integrand = y.*besselj(0,zeroes(i)*x).*x;  
    integral = h/3*( integrand(1) + 4*sum(integrand(2:2:M-1))...
        + 2*sum(integrand(3:2:M-2)) + integrand(M) ) ;
    Ck(i) = 2/besselj(1,zeroes(i))^2*integral;
    
end

% Display the hidden message

disp(decipher(Ck))