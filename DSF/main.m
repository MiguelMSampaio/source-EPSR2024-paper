% DSF modeling with Newton Rapshon power flow
% Author: Miguel Sampaio

clear
close all
clc

tInicio = clock;

%% ========================  DATA AND PRE PROCESSING =================================
% ------- Constants
it_max = 50;
ps.tol = 1e-6;        
ps.Sbase = 100e6;      %100 MVA

it_maxCtrl = 50;
ps.tolCtrl = 1e-6;

% ------- Input data

mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1);
filename = "data30_PVSysS";            % File .m with test system input data
run(strcat(newdir,"\test_systems\",filename,".m"))     

% ------- Structs definition
defVar_PU

%% ======================== NEWTON-RAPHSON METHOD =============================
it = 1;
itCtrl = 1;
itAcc = 0;

passoV = zeros(fv.NinvCtrl,1);
passoVold = zeros(fv.NinvCtrl,1);
criQ = zeros(fv.NinvCtrl,1);
criP = zeros(fv.NinvCtrl,1);
erroPrint = [];
erroGeral = [];

t0 = clock;
aux_index = sort([DBAR.index_PV' DBAR.index_PQ']);

% --- State variables initialization
X = [zeros(Nb-1,1); ones(DBAR.NPQ,1)];

% --- PV system power injection calculation
Processo1              



%% ========== OUTER LOOP (CONTROL AND LIMITS)
while(itCtrl<it_maxCtrl)  
    % --- First mismatch calculation
    [Pcalc,Qcalc,soma1,soma3] = calcPot(DBAR,DRAM,Nb,Nramo);       
    mismatch = [DBAR.Pesp;DBAR.Qesp] - [Pcalc;Qcalc];               
    aux = [ DBAR.index_ref ];                                       
    aux1 = [Nb+DBAR.index_ref', Nb+DBAR.index_PV'];                 
    mismatch([aux aux1]) = [];

    erro = max(abs(mismatch));                

    erroPrint(it) = erro;
    
    %% =========== INNER LOOP (POWER FLOW)
    while(erro>ps.tol && it<it_max)

        [J] = Jacobiano2(DBAR,DRAM,Nb,Nramo,Pcalc,Qcalc,soma1,soma3);

        deltaX = -J\mismatch;                                        
                                                                    
        X = X + deltaX;                                             

        DBAR.theta(aux_index) = X(1:(DBAR.NPQ+DBAR.NPV));
        DBAR.V(DBAR.index_PQ') = X(DBAR.NPQ+DBAR.NPV+1:2*DBAR.NPQ+DBAR.NPV);

        % mismatch
        [Pcalc,Qcalc,soma1,soma3] = calcPot(DBAR,DRAM,Nb,Nramo);        
        mismatch = [DBAR.Pesp;DBAR.Qesp] - [Pcalc;Qcalc];               
        aux = [ DBAR.index_ref ];                                       
        aux1 = [Nb+DBAR.index_ref', Nb+DBAR.index_PV'];                 
        mismatch([aux aux1]) = [];

        erro = max(abs(mismatch));                   

        it = it + 1;

        erroPrint(it) = erro;
    end
    %% ======== POWER FLOW END
    if it==it_max
        figure
        semilogy(1:it, erroPrint,'-x')
        title('Erro x Iteração')
        grid on
        error("Fluxo não convergiu em %d iterações, na %d itCtrl",it,itCtrl)
    end
    
    itAcc = itAcc + it;
    erroGeral(end+1:end+it) = erroPrint;
    it = 1;
    
    clear erroPrint
    
    %% ========= INVCONTROL
    if fv.NinvCtrl > 0
        
        % -- Get monitored voltage and P and Q criteria values
        for i = 1:fv.NinvCtrl
            fvI = fv.index_invCtrl(i);       
            barI = fv.bar(fvI);                
            
            fv.Vmon(i,itCtrl+1) = DBAR.V(barI);
%             fprintf('Vmon = %.7f\n',fv.Vmon(i,itCtrl+1))
            
            if itCtrl>1
                % Monitored voltage step
                passoV(i) = fv.Vmon(i,itCtrl+1) - fv.Vmon(i,itCtrl); 
                % Q criteria
                criQ(i) = abs(fv.Qd_end(i) - fv.Qca(fvI));
                % P criteria
                if fv.Plimit(fvI) > fv.Pl_end(i)
                    criP(i) = abs(fv.Plimit(fvI) - fv.Pl_end(i));
                else
                    criP(i) = 0;
                end
            end
        end
        
        % -- Convergence check
        if itCtrl>1
            if max(abs(passoV))<ps.tolCtrl && max(abs(criQ))<ps.tolCtrl && max(abs(criP))<ps.tolCtrl
%             if max(abs(passoV))<ps.tolCtrl
                break
            end
        end
        
        % -- Process 2 (control actions)
        Processo2
        
        passoVold = passoV;
        
        % -- Process 1 
        Processo1
        
        itCtrl = itCtrl + 1;
    else
        break
    end
    
end
%% ======= CONTROL LOOP END     


if itCtrl==it_maxCtrl    
    figure
    semilogy(1:itAcc, erroGeral,'-x')
    title('Erro x Iteração')
    grid on
    error("Laço de controle não convergiu na %d itCtrl",itCtrl)
end

tf=clock;
tempoNR = etime(tf,tInicio);                     

%% ================== POWER CALCULATION FOR EVERY BUS
Pfinal = DBAR.Pesp;                    
Qfinal = DBAR.Qesp;
for i=1:DBAR.NPV                                % PV bus type Q calculation
    n = DBAR.index_PV(i);    
    Qfinal(n) = Qcalc(n);
    DBAR.Qgen(n) = Qfinal(n) + DBAR.Qc(n);     
end

Qfinal(DBAR.index_ref) = Qcalc(DBAR.index_ref) + DBAR.Qc(DBAR.index_ref);     % P and Q for slack bus
Pfinal(DBAR.index_ref) = Pcalc(DBAR.index_ref) + DBAR.Pc(DBAR.index_ref);
DBAR.Pgen(DBAR.index_ref) = Pfinal(DBAR.index_ref);      
DBAR.Qgen(DBAR.index_ref) = Qfinal(DBAR.index_ref);


%% ============= FLOWS AND LOSSES CALCULATION
Yserie = DRAM.y./DRAM.a;          
susceptancia = zeros(Nramo,2);
for k=1:Nramo       
    susceptancia(k,1) = DRAM.y(k)*(1/DRAM.a(k))*(1/DRAM.a(k)-1) + 1i*DRAM.Bsh(k);     
    susceptancia(k,2) = (1-1/DRAM.a(k))*DRAM.y(k) + 1i*DRAM.Bsh(k);
end

Vcomplex = DBAR.V.*(cos(DBAR.theta) + 1i*sin(DBAR.theta));     

Ss = Vcomplex(DRAM.de).*conj((Vcomplex(DRAM.de) - Vcomplex(DRAM.para)).*Yserie + Vcomplex(DRAM.de).*susceptancia(:,1));     % From bus flows
Sr = Vcomplex(DRAM.para).*conj((Vcomplex(DRAM.para) - Vcomplex(DRAM.de)).*Yserie + Vcomplex(DRAM.para).*susceptancia(:,2)); % To bus flows
Pkm = real(Ss);   
Qkm = imag(Ss);     
Pmk = real(Sr);   
Qmk = imag(Sr);

%------------ losses
P_perdas = Pkm + Pmk;      
perdasPTot = sum(P_perdas); 
Q_perdas = Qkm + Qmk;       
perdasQTot = sum(Q_perdas);

tFinal = clock;


%% ======================= RESULTS ================================
puSI = ps.Sbase/1e6;    % in MVA
theta_grau = DBAR.theta*(180/pi);

% --- Save results in reuslts.mat file
ps.Nb = Nb;
ps.results.gerais = [ itAcc, erro, etime(tFinal,tInicio), tempoNR, perdasPTot*puSI, perdasQTot*puSI, itCtrl];
ps.results.barra = [ DBAR.bar   DBAR.V   theta_grau    DBAR.Pgen*puSI   DBAR.Qgen*puSI  DBAR.Pc*puSI  DBAR.Qc*puSI];
ps.results.ramo = [ DRAM.i  DRAM.de  DRAM.para   Pkm*puSI    Qkm*puSI     Pmk*puSI     Qmk*puSI    P_perdas*puSI   Q_perdas*puSI ];
ps.results.conv = erroGeral;
ps.results.fv = fv;

save('results.mat','ps');


PrintResults(ps);
PrintGraficos(ps);

