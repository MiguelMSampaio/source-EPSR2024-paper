%% ===== FUNC FOR CALCULATE THE VALUE OF PTcurve GIVEN T ========
% determine the piece-wise linear section of PTcurve and calculate the value

function y = PTcurve(curva,T)
% ENTRADA:      curva - 1x8 with the curve data T1,T2,T3,T4 P1,P2,P3,P4
%               T - temperature
% SAIDA:        y - Pnom in pu

T1 = curva(1,1); T2 = curva(1,2); T3 = curva(1,3); T4 = curva(1,4);
P1 = curva(1,5); P2 = curva(1,6); P3 = curva(1,7); P4 = curva(1,8);

if T < T1
    y = P1;
elseif T1 <= T && T < T2
    y = P1 + ((P2-P1)/(T2-T1))*(T - T1);
elseif T2 <= T && T < T3
    y = P2 + ((P3-P2)/(T3-T2))*(T - T2);
elseif T3 <= T && T < T4
    y = P3 + ((P4-P3)/(T4-T3))*(T - T3);
elseif T4 <= T
    y = P4;
end
    

end