%% ============== VOLT-VAR - look-up table =================

function Q = FunVV(Vmon,pontosVV,Qbase)
%   ENTRADA     - Vmon: monitored voltage
%               - pontosVV: vector with VV curve points
%               - Qbase: base reactive power 
%
%   SAÍDA       - Q: reactive power of VV curve for Vmon

VL = 0.5;
V1 = pontosVV(1);
V2 = pontosVV(2);
V3 = pontosVV(3);
V4 = pontosVV(4);
VH = 1.5;
Q1 = pontosVV(5)*Qbase;
Q2 = pontosVV(6)*Qbase;
Q3 = pontosVV(7)*Qbase;
Q4 = pontosVV(8)*Qbase;

if Vmon>VL && Vmon<=V1       
    Q = Q1;
elseif Vmon>V1 && Vmon<=V2
    Q = Q1 + ((Q2-Q1)/(V2-V1))*(Vmon - V1);            
elseif Vmon>V2 && Vmon<=V3
    Q = 0;                                                   
elseif Vmon>V3 && Vmon<=V4
    Q = Q3 + ((Q4-Q3)/(V4-V3))*(Vmon - V3);       
elseif Vmon>V4 && Vmon<=VH
    Q = Q4;
else
    Q = 0;
end

