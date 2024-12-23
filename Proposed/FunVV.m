% ============== FUNÇÃO VOLT-VAR - look-up table =================
% Calcula a potência reativa para a função de controle Volt-Var e a atualiza no objeto

function Q = FunVV(Vmon,pontosVV)
%   ENTRADA     - Vmon: tensão monitorada
%               - pontosVV: matriz com os pontos da curva de controle volt-var
%
%   SAÍDA       - Q: potencia reativa da curva VV para Vmon

% ---- Configuração dos valores de tensão e pot reativa da curva IVV
VL = 0.5;
V1 = pontosVV(1,2);
V2 = pontosVV(2,2);
V3 = pontosVV(3,2);
V4 = pontosVV(4,2);
VH = 1.5;
Qmax = pontosVV(1,1);
Qmin = pontosVV(4,1);


if Vmon>=VL && Vmon<=V1       
    Q = Qmax;
elseif Vmon>V1 && Vmon<=V2
    Q = (-Qmax/(V2-V1))*Vmon + (Qmax*V2)/(V2-V1);            %Equação da reta 1
elseif Vmon>V2 && Vmon<=V3
    Q = 0;                                                      %Deadband
elseif Vmon>V3 && Vmon<=V4
    Q = (Qmin/(V4-V3))*Vmon - (Qmin*V3)/(V4-V3);        %Equação da reta 2    
elseif Vmon>V4 && Vmon<=VH
    Q = Qmin;
else
    Q = 0;
end

