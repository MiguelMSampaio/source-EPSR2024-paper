%% =============  MISMATCHES ===========
% --- Grid
%   only .Pesp and .Qesp for buses without PV
%   -Smk calculation for buses with PV
% --- PV system
%   deltaG - general model
%   deltaPhi - controls. Complementarity/Smoothing

function [mismatch,DBAR,fv] = Mismatch(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps)
%% === constant
% kfp = tan(acos(fpLim)); 
if ~isempty(fv)
    xi = (ps.tol)^2;
    V1 = fv.VV(2);
    V2 = fv.VV(4);
    V3 = fv.VV(6);
    V4 = fv.VV(8);
end
    
%% =====  Skm, Smk, Qdisp and VV curve definition
DBAR.Pfv = zeros(Nb,1);
DBAR.Qfv = zeros(Nb,1);
if ~isempty(fv)
    if fv.Nvv > 0
        fv.splines(fv.Nvv,4) = SplineSuav;
        kVV = 1;
    end
end

for i=1:Nfv   
    barFV = fv.bar(i);  
    
    fv.Pkm(i) = (3/8)*fv.m(i)^2*fv.Vfv(i)^2*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
    fv.Pmk(i) = DBAR.V(barFV)^2*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
    
    fv.Qkm(i) = - (3/8)*fv.m(i)^2*fv.Vfv(i)^2*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
    fv.Qmk(i) = - DBAR.V(barFV)^2*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
    
    % sum, in case have more than one PV system on the same bus
    DBAR.Pfv(barFV) = DBAR.Pfv(barFV) + fv.Pkm(i);      %  do not consider transformer losses
    DBAR.Qfv(barFV) = DBAR.Qfv(barFV) + fv.Qkm(i);
%     DBAR.Pfv(barFV) = DBAR.Pfv(barFV) - fv.Pmk(i);      % with transformer losses
%     DBAR.Qfv(barFV) = DBAR.Qfv(barFV) - fv.Qmk(i);

    % - Available reactive power
    if abs(fv.Pkm(i))>=fv.Snom(i)
        fv.Qdisp(i,1) = 1e-1;     
        fprintf('Qdisp fixado\n')
    else
        fv.Qdisp(i,1) = sqrt(fv.Snom(i)^2 - fv.Pkm(i)^2);
    end    
%     fv.Qdisp(i,1) = sqrt(fv.Snom(i)^2 - fv.Pkm(i)^2);

    % --- define splines curves for smoothing for this it (Qdisp)
    if fv.cont(i) == 3
        Q1 = fv.VV(3)*fv.Qdisp(i,1);
        Q2 = fv.VV(5)*fv.Qdisp(i,1);
        Q3 = fv.VV(7)*fv.Qdisp(i,1);
        Q4 = fv.VV(9)*fv.Qdisp(i,1);

        fv.splines(kVV,1) = SplineSuav(V1-0.5,V1,V2, Q1,Q1,Q2, ps);
        fv.splines(kVV,2) = SplineSuav(V1,V2,V3, Q1,Q2,Q3, ps);
        fv.splines(kVV,3) = SplineSuav(V2,V3,V4, Q2,Q3,Q4, ps);
        fv.splines(kVV,4) = SplineSuav(V3,V4,V4+0.5, Q3,Q4,Q4, ps);
        
        kVV = kVV + 1;
    end

end


DBAR.Pesp = DBAR.Pgen - DBAR.Pc + DBAR.Pfv;
DBAR.Qesp = DBAR.Qgen - DBAR.Qc + DBAR.Qfv;


%%  =====  Grid mimsmatch calc

mismatch_grid = [DBAR.Pesp;DBAR.Qesp] - [Pcalc;Qcalc];      % Using Sesp in DBAR

aux = [ DBAR.index_ref ];                         
aux1 = [Nb+DBAR.index_ref', Nb+DBAR.index_PV'];   
mismatch_grid([aux aux1]) = [];


%% ==== PV system mismatch calc

if ~isempty(fv)
    % --- Every PV system
    deltaG1 = zeros(Nfv,1);   deltaG2 = zeros(Nfv,1); 
    deltaG3 = zeros(Nfv,1); 
    % Inverter capacity limit
    deltaPhi1 = zeros(Nfv,1);   deltaPhi2 = zeros(Nfv,1);
    
    % --- FPF
    deltaGfpf = zeros(fv.Nfpf,1);  
    deltaPhiFPF1 = zeros(fv.Nfpf,1);
    deltaPhiFPF2 = zeros(fv.Nfpf,1); deltaPhiFPF3 = zeros(fv.Nfpf,1);
    kFPF = 1;  
    
    % --- VW
    h = zeros(fv.Nvw,1);
    deltaPhiVW1 = zeros(fv.Nvw,1);  deltaPhiVW2 = zeros(fv.Nvw,1); 
    deltaPhiVW3 = zeros(fv.Nvw,1);  deltaPhiVW4 = zeros(fv.Nvw,1);  
    deltaPhiVW5 = zeros(fv.Nvw,1);  deltaPhiVW6 = zeros(fv.Nvw,1);      
    kVW = 1;
    
    % --- VV
    deltaPhiVV = zeros(fv.Nvv,1);  kVV = 1;
    
    for i=1:Nfv        
        barFV = fv.bar(i);    
    
        % base
        deltaG1(i) = fv.Iph(i) - fv.Is(i)*(exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) - 1) - fv.Ifv(i);
        deltaG2(i) = fv.Ifv(i) - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i)))/(fv.a(i)*fv.Vt(i));
        deltaG3(i) = fv.eta(i)*fv.Vfv(i)*fv.Ifv(i) - fv.Pkm(i);
        
        % inverter limit
        deltaPhi1(i) = deltaG2(i) + fv.xiC(i);
        deltaPhi2(i) = FunFBrelax((fv.Snom(i)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.xiC(i), xi);
        

        if fv.cont(i) == 1      % FPF
            deltaGfpf(kFPF) = fv.Qkm(i) - fv.Pkm(i)*tan(acos(fv.FPF(2)));
            
            deltaPhiFPF1(kFPF) = deltaGfpf(kFPF) - fv.Qli(kFPF) + fv.Qls(kFPF);
            deltaPhiFPF2(kFPF) = FunFBrelax((fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi);
            deltaPhiFPF3(kFPF) = FunFBrelax((fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi);
            
            kFPF = kFPF + 1;
        elseif fv.cont(i) == 2  % VW
            P1 = fv.VW(3)*fv.Snom(i); P2 = fv.VW(5)*fv.Snom(i);
            V1 = fv.VW(2);  V2 = fv.VW(4);
            h(kVW) = P1 + ((P2 - P1)/(V2 - V1))*(DBAR.V(barFV) - V1); 
            
            deltaPhiVW1(kVW) = deltaPhi1(i) + fv.Pls(kVW);
            deltaPhiVW2(kVW) = fv.Pmax(kVW) - h(kVW) + fv.PmaxLS(kVW) - fv.PmaxLI(kVW);
            
            deltaPhiVW3(kVW) = FunFBrelax((fv.Pmax(kVW)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.Pls(kVW), xi);
            deltaPhiVW4(kVW) = FunFBrelax((P1-fv.Pmax(kVW)), fv.PmaxLS(kVW), xi);
            deltaPhiVW5(kVW) = FunFBrelax((fv.Pmax(kVW)-P2), fv.PmaxLI(kVW), xi);
            
            deltaPhiVW6(kVW) = fv.Qkm(i);
            
            kVW = kVW + 1;
        elseif fv.cont(i) == 3          % VV 
            Q1 = fv.VV(3)*fv.Qdisp(i);
            Q2 = fv.VV(5)*fv.Qdisp(i);
            Q3 = fv.VV(7)*fv.Qdisp(i);
            Q4 = fv.VV(9)*fv.Qdisp(i);
            Vm = DBAR.V(fv.bar(i));
            
            if Vm <= fv.splines(kVV,1).v1
                deltaPhiVV(kVV) = fv.Qkm(i) - Q1;
            elseif fv.splines(kVV,1).v1 < Vm && Vm <= fv.splines(kVV,1).v2
                Qs = output(fv.splines(kVV,1),Vm);
                deltaPhiVV(kVV) = fv.Qkm(i) - Qs;
            elseif fv.splines(kVV,1).v2 < Vm && Vm <= fv.splines(kVV,2).v1
                deltaPhiVV(kVV) = fv.Qkm(i) - (Q1 + ((Q2-Q1)/(V2-V1))*(Vm-V1));
            elseif fv.splines(kVV,2).v1 < Vm && Vm <= fv.splines(kVV,2).v2
                Qs = output(fv.splines(kVV,2),Vm);
                deltaPhiVV(kVV) = fv.Qkm(i) - Qs;
            elseif fv.splines(kVV,2).v2 < Vm && Vm <= fv.splines(kVV,3).v1
                deltaPhiVV(kVV) = fv.Qkm(i);
            elseif fv.splines(kVV,3).v1 < Vm && Vm <= fv.splines(kVV,3).v2
                Qs = output(fv.splines(kVV,3),Vm);
                deltaPhiVV(kVV) = fv.Qkm(i) - Qs;
            elseif fv.splines(kVV,3).v2 < Vm && Vm <= fv.splines(kVV,4).v1
                deltaPhiVV(kVV) = fv.Qkm(i) - (Q4 + ((Q4-Q3)/(V4-V3))*(Vm-V4));
            elseif fv.splines(kVV,4).v1 < Vm && Vm <= fv.splines(kVV,4).v2
                Qs = output(fv.splines(kVV,4),Vm);
                deltaPhiVV(kVV) = fv.Qkm(i) - Qs;
            elseif fv.splines(kVV,4).v2 < Vm 
                deltaPhiVV(kVV) = fv.Qkm(i) - Q4;
            end
            
            kVV = kVV + 1;
        end
    end

    mismatch_fv = [ deltaG1; deltaG3; ...
        deltaPhi1(sort([fv.index_FPF' fv.index_VV'])); deltaPhi2; ...
        deltaPhiFPF1; deltaPhiFPF2; deltaPhiFPF3; ...
        deltaPhiVW1; deltaPhiVW2; deltaPhiVW3; deltaPhiVW4; deltaPhiVW5; deltaPhiVW6; ...
        deltaPhiVV ];

    mismatch = [mismatch_grid; mismatch_fv];
else
    mismatch = mismatch_grid;
end

end