% ============== FUNÇÃO FISCHER-BURMENSTEIN RELAXADA =================
% Calcula o resultado da funcao FB relaxda 

function Y = FunFBrelax(a,b,xi)
%   ENTRADA     - a,b : parametros da funcao FB
%               - xi: elemento de relaxamento da funcao FB
%
%   SAÍDA       - Y : resultado da funcao FB relaxada

Y = sqrt(a^2 + b^2 + xi) - (a + b);

end

