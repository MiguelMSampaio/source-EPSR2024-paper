% ============== VOLT-WATT - look-up table =================

function Plim = FunVW(Vmon,pontosVW,Pbase)
%   ENTRADA     - Vmon: monitored voltage
%               - pontosVW: vector with VW curve points
%               - Pbase: base reactive power
%
%   SAÍDA       - Plim: active power limit of VW curve for Vmon

VL = 0.5;
V1 = pontosVW(1);
V2 = pontosVW(2);
VH = 2;
P1 = pontosVW(3)*Pbase;
P2 = pontosVW(4)*Pbase;

if VL<Vmon && Vmon<=V1       
    Plim = P1;
elseif V1<Vmon && Vmon<=V2
    Plim = P1 + ((P2-P1)/(V2-V1))*(Vmon-V1);         
elseif V2<Vmon && Vmon<=VH
    Plim = P2;
else
    Plim = 0;
end

end

