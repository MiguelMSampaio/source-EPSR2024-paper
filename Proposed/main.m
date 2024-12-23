% Proposed modeling with Newton Raphson power flow
% Author: Miguel Sampaio

clear
close all
clc

tInicio = clock;

%% ========================  DATA =================================
% ----- Constants
it_max = 50;               % max nº it
ps.tol = 1e-6;              % tolerance mismatches in pu

it_parc = 6;               % nº it to partial result
tol_parc = 1e-2;            % tolerance to partial result

k_boltz = 1.3806503e-23;    % constant boltzmann
q = 1.60217646e-19;         % electron charge
Tn = 25; Gn = 1000;         % STC

% ----- Input data
mydir  = pwd;
idcs   = strfind(mydir,'/');
newdir = mydir(1:idcs(end)-1);
filename = 'data30_FVs';            
run(strcat(newdir,"\test_systems\",filename,".m"))     

% ---- Structs definition
defVar_PU

% ---- Func to calculate the 3 parameters (ideal model)
if ~isempty(fv)
    
    % 3 parameters in SI at module level
    for i = 1:Nfv
        fv.Iph(i) = (fv.Isc(i) + fv.Ki(i)*(fv.Tin(i) - Tn))*(fv.Gin(i)/Gn);
        fv.Vt(i) = fv.ns(i)*k_boltz*(fv.Tin(i)+273.15)/q;
        fv.Is(i) = ( fv.Isc(i) + fv.Ki(i)*(fv.Tin(i) - Tn) ) / ( exp((fv.Voc(i) + fv.Kv(i)*(fv.Tin(i)-Tn))/(fv.a(i)*fv.Vt(i))) );
    end

    % PV array
    fv.Iph = fv.Iph .* fv.Npp;
    fv.Vt = fv.Vt .* fv.Nss;
    fv.Is = fv.Is .* fv.Npp;

    % get in PU
    for i=1:Nfv
        fv.Iph(i) = fv.Iph(i) / (ps.Sbase/fv.Vbase(i));
        fv.Vt(i) = fv.Vt(i) / fv.Vbase(i);
        fv.Is(i) = fv.Is(i) / (ps.Sbase/fv.Vbase(i));
    end
end

%% ===========  State variables initialization
% electric grid (flat start)
Xg = [zeros(Nb-1,1); ones(DBAR.NPQ+DBAR.NFV,1)];       
% via input data file
% Xg = [zeros(Nb-1,1); DBAR.V(sort([DBAR.index_PQ' DBAR.index_FV']))]; 

% --- internal variable PV
if ~isempty(fv)
    % - Via Voc and Isc
%     fv.Vfv = fv.Voc .* fv.Nss ./ fv.Vbase;                  
%     fv.Ifv = fv.Isc .* fv.Npp ./ (ps.Sbase./fv.Vbase); 
    % - determinied values
%     fv.Vfv = ones(Nfv,1);
%     fv.Ifv = zeros(Nfv,1);
    % - Search for near MPP
    [fv.Vfv,fv.Ifv] = InicializaMPP(fv,Nfv);

    fv.m = ones(Nfv,1);                % ma unitary
    fv.alpha = zeros(Nfv,1);            % alphas zeros
    
    % slack variables and Pmax of VW
    fv.xiC = zeros(Nfv,1);                 
    fv.Qli = zeros(fv.Nfpf,1);             
    fv.Qls = zeros(fv.Nfpf,1);             
    fv.Pmax = fv.Snom(fv.index_VW,1).*ones(fv.Nvw,1); 
    fv.Pls = zeros(fv.Nvw,1);              
    fv.PmaxLI = zeros(fv.Nvw,1);           
    fv.PmaxLS = zeros(fv.Nvw,1);          
   
    Xpv = [fv.Vfv; fv.Ifv; fv.m; fv.alpha];
    Xfolga = [fv.xiC; fv.Qli; fv.Qls; fv.Pmax; fv.Pls; fv.PmaxLI; fv.PmaxLS];
    
    X = [Xg; Xpv];
else
    X = Xg;
end

%% ======================== NEWTON-RAPHSON METHOD =============================
it = 0;
t0 = clock;
flag_div = 0;

aux_index = sort([DBAR.index_PV' DBAR.index_PQ' DBAR.index_FV']);

% ---- First mismatch calc
[Pcalc,Qcalc,soma1,soma3] = CalcPot(DBAR,DRAM,Nb,Nramo);    
[mismatch,DBAR,fv] = MismatchIni2(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps); 
erro = max(abs(mismatch));    

erroPrint(it+1) = erro;

% Data to convergence analysis of internal variables, only of first PV system
ColetaDadosConvIni

% ---- System 1
while(erro>tol_parc && it<it_parc)
    
    [J] = JacobianoIni2(DBAR,DRAM,fv,Nb,Nramo,Nfv,Pcalc,Qcalc,soma1,soma3);
    
    deltaX = -J\mismatch;               
    
    X = X + deltaX;     
    
    DBAR.theta(aux_index) = X(1:(DBAR.NPQ+DBAR.NPV+DBAR.NFV));
    DBAR.V(sort([DBAR.index_PQ' DBAR.index_FV'])) = X(DBAR.NPQ+DBAR.NPV+DBAR.NFV+1:2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV);
    if ~isempty(fv)
        aux_end = 2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV;
        fv.Vfv = X(aux_end + 1:aux_end+Nfv);
        fv.Vfv = X((aux_end + 1):(aux_end + Nfv));
        fv.Ifv = X((aux_end + Nfv+1):(aux_end + 2*Nfv));
        fv.m = X((aux_end + 2*Nfv+1):(aux_end + 3*Nfv));
        fv.alpha = X((aux_end + 3*Nfv+1):(aux_end + 4*Nfv));
    end
    
    [Pcalc,Qcalc,soma1,soma3] = CalcPot(DBAR,DRAM,Nb,Nramo);
    [mismatch,DBAR,fv] = MismatchIni2(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps);
    
    erro = max(abs(mismatch));      

    it = it + 1;
    erroPrint(it+1) = erro;
    
    ColetaDadosConvIni
end

it_resParc = it;
if ~isempty(fv)
    X = [X; Xfolga];
end

% mismatch for the new state variables
[Pcalc,Qcalc,soma1,soma3] = CalcPot(DBAR,DRAM,Nb,Nramo);
[mismatch,DBAR,fv] = Mismatch(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps);

erro = max(abs(mismatch));  

% ---- System 2
while(erro>ps.tol && it<it_max)

    [J] = Jacobiano(DBAR,DRAM,fv,Nb,Nramo,Nfv,Pcalc,Qcalc,soma1,soma3,ps);      
    
    deltaX = -J\mismatch;     
    
    X = X + deltaX;           
    
    DBAR.theta(aux_index) = X(1:(DBAR.NPQ+DBAR.NPV+DBAR.NFV));            
    DBAR.V(sort([DBAR.index_PQ' DBAR.index_FV'])) = X(DBAR.NPQ+DBAR.NPV+DBAR.NFV+1:2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV);  
    if ~isempty(fv)
        aux_end = 2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV;
        fv.Vfv = X((aux_end + 1):(aux_end + Nfv));
        fv.Ifv = X((aux_end + Nfv+1):(aux_end + 2*Nfv));
        fv.m = X((aux_end + 2*Nfv+1):(aux_end + 3*Nfv));
        fv.alpha = X((aux_end + 3*Nfv+1):(aux_end + 4*Nfv));
        fv.xiC = X((aux_end + 4*Nfv+1):(aux_end + 5*Nfv));
        if fv.Nfpf>0
            fv.Qli = X((aux_end + 5*Nfv + 1):(aux_end + 5*Nfv + fv.Nfpf));
            fv.Qls = X((aux_end + 5*Nfv + fv.Nfpf + 1):(aux_end + 5*Nfv + 2*fv.Nfpf));
        end
        if fv.Nvw>0
            fv.Pmax = X((aux_end + 5*Nfv + 2*fv.Nfpf + 1):(aux_end+5*Nfv+2*fv.Nfpf + fv.Nvw));
            fv.Pls = X((aux_end + 5*Nfv + 2*fv.Nfpf + fv.Nvw + 1):(aux_end+5*Nfv+2*fv.Nfpf + 2*fv.Nvw));
            fv.PmaxLI = X((aux_end + 5*Nfv + 2*fv.Nfpf + 2*fv.Nvw + 1):(aux_end+5*Nfv+2*fv.Nfpf + 3*fv.Nvw));
            fv.PmaxLS = X((aux_end + 5*Nfv + 2*fv.Nfpf + 3*fv.Nvw + 1):(aux_end+5*Nfv+2*fv.Nfpf + 4*fv.Nvw));
        end
    end
        
    [Pcalc,Qcalc,soma1,soma3] = CalcPot(DBAR,DRAM,Nb,Nramo);
    [mismatch,DBAR,fv] = Mismatch(Pcalc,Qcalc,DBAR,fv,Nb,Nfv,ps);
    
    erro = max(abs(mismatch));  

    it = it + 1;
        
    erroPrint(it+1) = erro;
    
    ColetaDadosConv
    
     if max(abs(X))>1000
        flag_div = 1;
        break
    end
    
end


if ~isempty(fv) && Nfv==1
    fileID = fopen('debugTrace.txt','w');
    if fv.cont(1)==1
        fprintf(fileID,'%3s %10s %10s %10s %10s %13s %10s %10s %10s %10s %10s %10s %10s\n','it','Vfv(pu)','Ifv(pu)','Pfv(pu)','m',...
            'alpha(deg)','Vk(pu)','Pkm(pu)','Qkm(pu)','Qdisp','xiC','Qli','Qls');
        fprintf(fileID,'%3d %10.4f %10.4f %10.4f %10.4f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',dataDebug');
    elseif fv.cont(1)==2
        fprintf(fileID,'%3s %10s %10s %10s %10s %13s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n','it','Vfv(pu)','Ifv(pu)','Pfv(pu)','ma',...
            'alpha(deg)','Vk(pu)','Pkm(pu)','Qkm(pu)','xiC','Pmax','Pls','PmaxLI','PmaxLS','Vm');
        fprintf(fileID,'%3d %10.4f %10.4f %10.4f %10.4f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',dataDebug');
    elseif fv.cont(1)==3
        fprintf(fileID,'%3s %10s %10s %10s %10s %13s %10s %10s %10s %10s %10s %10s\n','it','Vfv(pu)','Ifv(pu)','Pfv(pu)','m',...
            'alpha(deg)','Vk(pu)','Pkm(pu)','Qkm(pu)','Qdisp','xiC','Vm');
        fprintf(fileID,'%3d %10.4f %10.4f %10.4f %10.4f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',dataDebug');
    end
    fclose(fileID);
end

if it==it_max || flag_div==1
    figure
    semilogy(0:it, erroPrint,'-x')
    title('Erro x Iteração')
    grid on
    if ~isempty(fv) && Nfv==1
        figure 
        cmap = colormap(turbo(size(resPrint,2)));
        for k=1:size(resPrint,2)
            semilogy(0:it, abs(resPrint(:,k)),'Color', cmap(k, :))
            hold on
        end
        hold off
        legend('deltaG1','deltaG3')
        if fv.cont(1)==1
            legend('deltaG1','deltaG3','deltaPhi1','deltaPhi2','deltaPhiFPF1','deltaPhiFPF2','deltaPhiFPF3')
        elseif fv.cont(1)==2
            legend('deltaG1','deltaG3','deltaPhi2','deltaPhiVW1','deltaPhiVW2','deltaPhiVW3','deltaPhiVW4','deltaPhiVW5','deltaPhiVW6')
        elseif fv.cont(1)==3
            legend('deltaG1','deltaG3','deltaPhi1','deltaPhi2','deltaPhiVV')
        end
        
        figure
        cmap = colormap(turbo(size(varPrint,2)));
        hold on
        for k=1:size(varPrint,2)
            plot(0:it, varPrint(:,k),'Color',cmap(k,:))
        end
        hold off
        legend('Vfv','Ifv','m','alpha')
        if fv.cont(1)==1
            legend('Vfv','Ifv','m','alpha','xiC','Qli','Qls')
        elseif fv.cont(1)==2
            legend('Vfv','Ifv','m','alpha','xiC','Pmax','Pls','PmaxLI','PmaxLS')
        elseif fv.cont(1)==3
            legend('Vfv','Ifv','m','alpha','xiC')
        end
    end
    
    if it==it_max && flag_div==0
        error("Não convergiu no nº de iterações máximo")
    elseif flag_div == 1
        error("Divergiu")
    end
end

tf=clock;
tempoNR = etime(tf,tInicio);         



%% =======================  POWER CALCULATION FOR EVERY BUS ====================================
Pfinal = DBAR.Pesp;            
Qfinal = DBAR.Qesp;
for i=1:DBAR.NPV                      
    n = DBAR.index_PV(i);    
    Qfinal(n) = Qcalc(n);
    DBAR.Qgen(n) = Qfinal(n) + DBAR.Qc(n); 
end

Qfinal(DBAR.index_ref) = Qcalc(DBAR.index_ref) + DBAR.Qc(DBAR.index_ref);  
Pfinal(DBAR.index_ref) = Pcalc(DBAR.index_ref) + DBAR.Pc(DBAR.index_ref);
DBAR.Pgen(DBAR.index_ref) = Pfinal(DBAR.index_ref);  
DBAR.Qgen(DBAR.index_ref) = Qfinal(DBAR.index_ref);


%% ======================== FLOWS AND LOSSES CALCULATION ===================
Yserie = DRAM.y./DRAM.a;        
susceptancia = zeros(Nramo,2);
for k=1:Nramo       
    susceptancia(k,1) = DRAM.y(k)*(1/DRAM.a(k))*(1/DRAM.a(k)-1) + 1i*DRAM.Bsh(k);  
    susceptancia(k,2) = (1-1/DRAM.a(k))*DRAM.y(k) + 1i*DRAM.Bsh(k);
end

Vcomplex = DBAR.V.*(cos(DBAR.theta) + 1i*sin(DBAR.theta)); 

Ss = Vcomplex(DRAM.de).*conj((Vcomplex(DRAM.de) - Vcomplex(DRAM.para)).*Yserie + Vcomplex(DRAM.de).*susceptancia(:,1));   
Sr = Vcomplex(DRAM.para).*conj((Vcomplex(DRAM.para) - Vcomplex(DRAM.de)).*Yserie + Vcomplex(DRAM.para).*susceptancia(:,2)); 
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
for i=1:Nfv
    fv.fp(i,1) = sign(fv.Qkm(i)) * cos(atan(fv.Qkm(i)/fv.Pkm(i)));
end
if ~isempty(fv) && Nfv==1
    if ~isreal(dataDebug)
        fprintf('\n=== ! PROBLEMA NÚMERO COMPLEXO !\n')
    end
end

% --- Save results in reuslts.mat file
ps.Nb = Nb;
ps.fv = fv;
ps.results.gerais = [ it, erro, etime(tFinal,tInicio), tempoNR, perdasPTot*puSI, perdasQTot*puSI, it_resParc];
ps.results.barra = [ DBAR.bar   DBAR.V   theta_grau    (DBAR.Pfv+DBAR.Pgen)*puSI   (DBAR.Qfv+DBAR.Qgen)*puSI  DBAR.Pc*puSI  DBAR.Qc*puSI];
ps.results.ramo = [ DRAM.i  DRAM.de  DRAM.para   Pkm*puSI    Qkm*puSI     Pmk*puSI     Qmk*puSI    P_perdas*puSI   Q_perdas*puSI ];
ps.results.fv.Nfv = Nfv;
if Nfv>0
    ps.results.fv.varInt = [ fv.bar fv.Vfv fv.Ifv fv.Vfv.*fv.Ifv fv.m fv.alpha*(180/pi) ...
            sqrt(3/8)*fv.m.*fv.Vfv fv.Pkm fv.Qkm sqrt(fv.Pkm.^2+fv.Qkm.^2) fv.fp];
    ps.results.fv.varFolgaFPF = [fv.bar(fv.index_FPF,1) fv.xiC(fv.index_FPF,1) fv.Qli fv.Qls];
    ps.results.fv.varFolgaVW = [fv.bar(fv.index_VW,1) fv.xiC(fv.index_VW,1) fv.Pmax fv.Pls fv.PmaxLI fv.PmaxLS];
    ps.results.fv.varFolgaVV = [fv.bar(fv.index_VV,1) fv.xiC(fv.index_VV,1) ];
end
ps.results.conv = erroPrint;
if ~isempty(fv) && Nfv==1
    ps.results.resPrint = resPrint;
    ps.results.varPrint = varPrint;
end

save('results.mat','ps');

PrintResults(ps);
PrintGraficos(ps);
