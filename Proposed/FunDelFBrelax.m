% ============== DERIVADAS DA FUNÇÃO FISCHER-BURMENSTEIN RELAXADA =================
% Calcula o resultado das duas derivadas posiveis da funcao FB relaxda 
% formato de FB = sqrt(a(x)^2 + b^2 + xi) - (a + b)

function Y = FunDelFBrelax(da,a,b,xi,flag)
%   ENTRADA     - da: derivada de a 
%               - a,b : parametros da funcao FB
%               - xi: elemento de relaxamento da funcao FB
%               - flag: indica se eh derivada em relacao a 'x' (1) ou 'b' (2)
%
%   SAÍDA       - Y : resultado da derivada

if flag == 1
    Y = da*( (a/sqrt(a^2 + b^2 + xi)) - 1 );
elseif flag ==2
    Y = ( b/(sqrt(a^2 + b^2 + xi)) ) - 1;
end

end

