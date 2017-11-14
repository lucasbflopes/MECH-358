function [ code_deciphered ] = decipher( input )

% This code converts the array input filled with real numbers into the string 
% array code_deciphered. The rule of translating is 0->'space', 1->a, 2->b
% ,... and so on. Since "a" is 97 in ASCII, in order to fulfill this rule,
% we have to sum 96 to each component of input.

code = round(input)+96;

for i = 1:length(code)
    if code(i)==96
        code(i) = 32; % 32 is SPACE in ASCII
    end
end

code_deciphered = char(code);

end

