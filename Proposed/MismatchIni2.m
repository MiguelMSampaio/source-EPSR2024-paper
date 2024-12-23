%% ============= MISMATCHES CALC ===========
% --- Grid
%   only .Pesp and .Qesp for buses without PV
%   -Smk calculation for buses with PV
% --- PV system
%   deltaG - general model
%   deltaPhi - controls. Complementarity/Smoothing

function [mismatch,DBAR,fv] = MismatchIni2(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps)
%% ===== Skm, Smk, Qdisp and VV curve definition

DBAR.Pfv = zeros(Nb,1);
DBAR.Qfv = zeros(Nb,1);

for i=1:Nfv      
    barFV = fv.bar(i);     
    % Skm
    fv.Pkm(i) = (3/8)*fv.m(i)^2*fv.Vfv(i)^2*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
%     fv.Pmk(i) = DBAR.V(barFV)^2*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
%         * ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
    
    fv.Qkm(i) = - (3/8)*fv.m(i)^2*fv.Vfv(i)^2*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
        * ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
%     fv.Qmk(i) = - DBAR.V(barFV)^2*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV)...
%         * ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
    
    % sum, in case have more than one PV system on the same bus
    DBAR.Pfv(barFV) = DBAR.Pfv(barFV) + fv.Pkm(i);      % do not consider transformer losses
    DBAR.Qfv(barFV) = DBAR.Qfv(barFV) + fv.Qkm(i);
%     DBAR.Pfv(barFV) = DBAR.Pfv(barFV) - fv.Pmk(i);      % with transformer losses
%     DBAR.Qfv(barFV) = DBAR.Qfv(barFV) - fv.Qmk(i);

    % - Available reactive power
    if fv.Pkm(i)>=fv.Snom(i)
        fv.Qdisp(i,1) = 1e-1;       
        fprintf('Qdisp fixado\n')
    else
        fv.Qdisp(i,1) = sqrt(fv.Snom(i)^2 - fv.Pkm(i)^2);
    end    

end


DBAR.Pesp = DBAR.Pgen - DBAR.Pc + DBAR.Pfv;
DBAR.Qesp = DBAR.Qgen - DBAR.Qc + DBAR.Qfv;


%%  ===== Grid mimsmatch calc

mismatch_grid = [DBAR.Pesp;DBAR.Qesp] - [Pcalc;Qcalc];      % Using Sesp DBAR

aux = [ DBAR.index_ref ];                              
aux1 = [Nb+DBAR.index_ref', Nb+DBAR.index_PV'];        
mismatch_grid([aux aux1]) = [];


%% ==== PV system mismatch calc

if ~isempty(fv)
    % --- Every PV system
    deltaG1 = zeros(Nfv,1);   deltaG2 = zeros(Nfv,1); 
    deltaG3 = zeros(Nfv,1); 
    
    % --- FPF control
    deltaGfpf = zeros(fv.Nfpf,1);  
    kFPF = 1;
    
    % --- VW control
    deltaGvw = zeros(fv.Nvw,1);      
    kVW = 1;
    
    % --- VV control
    deltaGvv = zeros(fv.Nvv,1);  
    kVV = 1;
    
    for i=1:Nfv        
%         barFV = fv.bar(i);      % barra do FV
    
        % base
        deltaG1(i) = fv.Iph(i) - fv.Is(i)*(exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) - 1) - fv.Ifv(i);
        deltaG2(i) = fv.Ifv(i) - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i)))/(fv.a(i)*fv.Vt(i));
        deltaG3(i) = fv.eta(i)*fv.Vfv(i)*fv.Ifv(i) - fv.Pkm(i);

        % control
        if fv.cont(i) == 1      %  FPF
            deltaGfpf(kFPF) = fv.Qkm(i) - fv.Pkm(i)*tan(acos(fv.FPF(2)));
            
            kFPF = kFPF + 1;
        elseif fv.cont(i) == 2  %  VW
            deltaGvw(kVW) = fv.Qkm(i);
            
            kVW = kVW + 1;
        elseif fv.cont(i) == 3 %  VV (%Snom and priority watt)
            deltaGvv(kVV) = fv.Qkm(i);
            
            kVV = kVV + 1;
        end
    end

    mismatch_fv = [ deltaG1; deltaG2; deltaG3; ...
        deltaGfpf; deltaGvw; deltaGvv];

    mismatch = [mismatch_grid; mismatch_fv];
else
    mismatch = mismatch_grid;
end

end