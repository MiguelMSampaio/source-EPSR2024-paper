%% ================== JACOBIAN ===========================

function [J] = Jacobiano(DBAR,DRAM,fv,Nb,Nramo,Nfv,Pcalc,Qcalc,som1,som3,ps)
%% === constants
a = ones(Nramo,1)./DRAM.a;   
% kfp = tan(acos(fpLim));


%% == deltaSkm
if ~isempty(fv)
    xi = ps.tol^2;
    V1 = fv.VV(2);
    V2 = fv.VV(4);
    V3 = fv.VV(6);
    V4 = fv.VV(8);
    
    dPkm_dTm = zeros(Nfv,1);  dPkm_dVm = zeros(Nfv,1);
    dPkm_dVfv = zeros(Nfv,1);  
    dPkm_dm = zeros(Nfv,1);   
    dPkm_dalpha = zeros(Nfv,1);
    
    dQkm_dTm = zeros(Nfv,1);  dQkm_dVm = zeros(Nfv,1);
    dQkm_dVfv = zeros(Nfv,1);
    dQkm_dm = zeros(Nfv,1); dQkm_dalpha = zeros(Nfv,1);
    
    dPmk_dTm = zeros(Nfv,1); dPmk_dVm = zeros(Nfv,1);
    dPmk_dVfv = zeros(Nfv,1); dPmk_dm = zeros(Nfv,1);
    dPmk_dalpha = zeros(Nfv,1);
    dQmk_dTm = zeros(Nfv,1); dQmk_dVm = zeros(Nfv,1);
    dQmk_dVfv = zeros(Nfv,1); dQmk_dm = zeros(Nfv,1);
    dQmk_dalpha = zeros(Nfv,1);
    
    dQdisp_dTm = zeros(Nfv,1); dQdisp_dVm = zeros(Nfv,1);
    dQdisp_dVfv = zeros(Nfv,1);
    dQdisp_dm = zeros(Nfv,1); dQdisp_dalpha = zeros(Nfv,1);

    for i=1:Nfv
        barFV = fv.bar(i);   

        % -- dPkm
        dPkm_dTm(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) - fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
        dPkm_dVm(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
            ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
        dPkm_dVfv(i) = (3/4)*fv.m(i)^2*fv.Vfv(i)*fv.g_t(i) - sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
        dPkm_dm(i) = (3/4)*fv.m(i)*fv.Vfv(i)^2*fv.g_t(i) - sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
        dPkm_dalpha(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );    
        
        % -- dQkm
        dQkm_dTm(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );
        dQkm_dVm(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
            ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
        dQkm_dVfv(i) = - (3/4)*fv.m(i)^2*fv.Vfv(i)*fv.b_t(i) + sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
        dQkm_dm(i) = - (3/4)*fv.m(i)*fv.Vfv(i)^2*fv.b_t(i) + sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) + fv.b_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) );
        dQkm_dalpha(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*cos(fv.alpha(i) - DBAR.theta(barFV)) - fv.b_t(i)*sin(fv.alpha(i) - DBAR.theta(barFV)) );

        % -- dPmk
        dPmk_dTm(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
        dPmk_dVm(i) = 2*DBAR.V(barFV)*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
            ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        dPmk_dVfv(i) = - sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        dPmk_dm(i) = - sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        dPmk_dalpha(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) - fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
        
        % -- dQmk
        dQmk_dTm(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( - fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) - fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        dQmk_dVm(i) = -2*DBAR.V(barFV)*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
            ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
        dQmk_dVfv(i) = sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...                       
            ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
        dQmk_dm(i) = sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...                       
            ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
        dQmk_dalpha(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
            ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        
        % -- Qdisp
        dQdisp_dTm(i) = - dPkm_dTm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
        dQdisp_dVm(i) = - dPkm_dVm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
        dQdisp_dVfv(i) = - dPkm_dVfv(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
        dQdisp_dm(i) = - dPkm_dm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
        dQdisp_dalpha(i) = - dPkm_dalpha(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
            
    end
end

%% Grid

H = zeros(DBAR.NPQ+DBAR.NPV+DBAR.NFV, DBAR.NPQ+DBAR.NPV+DBAR.NFV);       % dDeltaP/dtheta
N = zeros(DBAR.NPQ+DBAR.NPV+DBAR.NFV, DBAR.NPQ+DBAR.NFV);       % dDeltaP/DV
M = zeros(DBAR.NPQ+DBAR.NFV, DBAR.NPQ+DBAR.NPV+DBAR.NFV);       % dDeltaQ/dtheta
L = zeros(DBAR.NPQ+DBAR.NFV, DBAR.NPQ+DBAR.NFV);       % dDeltaQ/dV

% == H - delDeltaP/delTheta (NPV+NPQ+NFV x NPV+NPQ+NFV) (Nb-1 x Nb-1)
for i=2:Nb
    for j=2:Nb
        [km,~] = km_flag(Nramo, DRAM, i, j);
        if i==j         
            H(i-1, j-1) = H(i-1, j-1) + Qcalc(i) + (DBAR.V(i)^2* (abs(DBAR.mhoSh(i)) + som3(i))); 
            for m = 1:Nfv   
                if j==fv.bar(m) 
                    H(i-1, j-1) = H(i-1, j-1) + dPkm_dTm(m); 
%                     H(i-1, j-1) = H(i-1, j-1) - dPmk_dTm(m); 
                end
            end
        else           
            if km~=0   
                H(i-1, j-1) = + DBAR.V(i)*DBAR.V(j)*a(km) * (DRAM.g(km)*sin(DBAR.theta(i) - DBAR.theta(j)) - ...
                    DRAM.b(km)*cos(DBAR.theta(i) - DBAR.theta(j)));                  
            end
        end
    end
end

auxIndex_PQ_FV = sort([DBAR.index_PQ' DBAR.index_FV']);

% == N - delDeltaP/delV (NPV+NPQ+NFV x NPQ+NFV)
for i=2:Nb
    for j=1:DBAR.NPQ+DBAR.NFV      
        n = auxIndex_PQ_FV(j);
        [km,~] = km_flag(Nramo, DRAM, i, n);
        
        if i==n    
            N(i-1, j) = N(i-1 ,j) - (Pcalc(i) + DBAR.V(i)^2*som1(i))/DBAR.V(i);
            for m = 1:Nfv   
                if n==fv.bar(m)  
                    N(i-1, j) = N(i-1, j) + dPkm_dVm(m); 
%                     N(i-1, j) = N(i-1, j) - dPmk_dVm(m);
                end
            end
        else       
            if km~=0                
                N(i-1, j) = + DBAR.V(i)*a(km) * (DRAM.g(km)*cos(DBAR.theta(i) - DBAR.theta(n)) + ...
                    DRAM.b(km)*sin(DBAR.theta(i) - DBAR.theta(n)));
            end
        end
    end
end


% == M  - delDeltaQ/delTheta (NPQ+NFV x NPV+NPQ+NFV)
for i=1:DBAR.NPQ+DBAR.NFV
    for j=2:Nb
        n = auxIndex_PQ_FV(i);
        [km,~] = km_flag(Nramo, DRAM, n, j);
        if n==j     
            M(i, j-1) = M(i, j-1) - Pcalc(n) + (DBAR.V(n)^2*som1(n));
            for m = 1:Nfv  
                if j==fv.bar(m)  
                    M(i, j-1) = M(i, j-1) + dQkm_dTm(m); 
%                     M(i, j-1) = M(i, j-1) - dQmk_dTm(m); 
                end
            end
        else      
            if km~=0
                M(i,j-1) = - DBAR.V(n)*DBAR.V(j)*a(km)*(DRAM.g(km)*cos(DBAR.theta(n) - DBAR.theta(j)) + ...
                    DRAM.b(km)*sin(DBAR.theta(n) - DBAR.theta(j)));
            end
        end
    end
end

% == L  - delDeltaQ/delV (NPQ+NFV x NPQ+NFV)
for i=1:DBAR.NPQ+DBAR.NFV
    n = auxIndex_PQ_FV(i);
    for j=1:DBAR.NPQ+DBAR.NFV
        k = auxIndex_PQ_FV(j);
        [km,~] = km_flag(Nramo, DRAM, n, k);
        if i==j             
            L(i,j) = L(i,j) - (Qcalc(n) - DBAR.V(n)^2*(abs(DBAR.mhoSh(n)) + som3(n))) / DBAR.V(n);
            for m = 1:Nfv   
                if k==fv.bar(m)
                    L(i,j) = L(i,j) + dQkm_dVm(m);  
%                     L(i,j) = L(i,j) - dQmk_dVm(m); 
                end
            end
        else                
            if km~=0
                L(i,j) = + DBAR.V(n)*a(km) * (DRAM.g(km)*sin(DBAR.theta(n) - DBAR.theta(k)) - ...
                    DRAM.b(km)*cos(DBAR.theta(n) - DBAR.theta(k)));
            end            
        end
    end
end

Jgrid = [H N; M L];


%% ==== Extensão do jacobiano

X = [];
Y = [];
Z = [];

if ~isempty(fv)    
    
    % X (2NPQ+2NFV+NPV  X  4Nfv+Nfv+2Nfv,fpf+4Nfv,vw) 
    %       delDeltaP e delDeltaQ / delVarInternas
    X = zeros(2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV, 5*Nfv+2*fv.Nfpf+4*fv.Nvw);

    for i=1:Nfv
        barFV = fv.bar(i);
        indexFV = find(auxIndex_PQ_FV==barFV);  
        

        % delDP/delVfv 
        X(barFV-1, i) = dPkm_dVfv(i);
        % delDP/delm
        X(barFV-1, 2*Nfv + i) = dPkm_dm(i);
        % delDP/delAlpha
        X(barFV-1, 3*Nfv + i) = dPkm_dalpha(i);

        % deltaDQ/deltaVfv
        X(DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV, i) = dQkm_dVfv(i);
        % deltaDQ/deltam
        X(DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV, 2*Nfv + i) = dQkm_dm(i);
        % deltaDQ/deltaAlpha
        X(DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV, 3*Nfv + i) = dQkm_dalpha(i);
    end

    
    % Y (4Nfv+Nfv+2Nfv,fpf+4Nfv,vw X 2NPQ+2Nfv+NPV) 
    %       delG/delTheta e delG/delVgrid
    Y = zeros(5*Nfv+2*fv.Nfpf+4*fv.Nvw, 2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV);
    
    kFPF=1; kVW=1; kVV=1;
    for i=1:Nfv
        barFV = fv.bar(i);
        indexFV = find(auxIndex_PQ_FV==barFV); 

        % deltaG3/Theta_m
        Y(Nfv + i, barFV-1) = -dPkm_dTm(i);
        % deltaG3/Vm
        Y(Nfv + i, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = -dPkm_dVm(i);
        
        
        if fv.cont(i) == 1          %  FPF      
            % deltaPhiFPF1/Theta_m
            Y(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, barFV-1) = dQkm_dTm(i) - dPkm_dTm(i)*tan(acos(fv.FPF(2)));
            % deltaPhiFPF1/Vm
            Y(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - dPkm_dVm(i)*tan(acos(fv.FPF(2)));
            
            % deltaPhiFPF2/Theta_m
            Y(3*Nfv+2*fv.Nfpf+fv.Nvv + kFPF, barFV-1) = FunDelFBrelax((dQkm_dTm(i)+dQdisp_dTm(i)), (fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi, 1);
            % deltaPhiFPF2/Vm
            Y(3*Nfv+2*fv.Nfpf+fv.Nvv + kFPF, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = FunDelFBrelax((dQkm_dVm(i)+dQdisp_dVm(i)), (fv.Qkm(i)+fv.Qdisp(i)),fv.Qli(kFPF), xi, 1);
            
            % deltaPhiFPF3/Theta_m
            Y(3*Nfv+3*fv.Nfpf+fv.Nvv + kFPF, barFV-1) = FunDelFBrelax((dQdisp_dTm(i)-dQkm_dTm(i)), (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi, 1);
            % deltaPhiFPF3/Vm
            Y(3*Nfv+3*fv.Nfpf+fv.Nvv + kFPF, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = FunDelFBrelax((dQdisp_dVm(i)-dQkm_dVm(i)), (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi,1);
            
            kFPF = kFPF+1;
            
        elseif fv.cont(i) == 2      %  VW
            P1 = fv.VW(3)*fv.Snom(i); P2 = fv.VW(5)*fv.Snom(i);
            V1 = fv.VW(2);  V2 = fv.VW(4);
            
            % deltaPhiVW2/Vm
            Y(3*Nfv+4*fv.Nfpf+fv.Nvv+fv.Nvw + kVW, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = -((P2 - P1)/(V2 - V1));
            
            % deltaPhiVW6/Theta_m
            Y(3*Nfv+4*fv.Nfpf+fv.Nvv+5*fv.Nvw + kVW, barFV-1) = dQkm_dTm(i);
            % deltaPhiVW6/Vm
            Y(3*Nfv+4*fv.Nfpf+fv.Nvv+5*fv.Nvw + kVW, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            
            kVW = kVW + 1;
            
        elseif fv.cont(i) == 3      %  VV
            Q1 = fv.VV(3)*fv.Qdisp(i);
            Q2 = fv.VV(5)*fv.Qdisp(i);
            Q3 = fv.VV(7)*fv.Qdisp(i);
            Q4 = fv.VV(9)*fv.Qdisp(i);
            Vm = DBAR.V(fv.bar(i));
            
            % deltaGvv/Theta_m
            Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, barFV-1) = dQkm_dTm(i);
            % deltaGvv/Vm
            if Vm <= fv.splines(kVV,1).v1
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            elseif fv.splines(kVV,1).v1 < Vm && Vm <= fv.splines(kVV,1).v2
                delQs = del(fv.splines(kVV,1),Vm);
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - delQs;
            elseif fv.splines(kVV,1).v2 < Vm && Vm <= fv.splines(kVV,2).v1
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - ((Q2-Q1)/(V2-V1));
            elseif fv.splines(kVV,2).v1 < Vm && Vm <= fv.splines(kVV,2).v2
                delQs = del(fv.splines(kVV,2),Vm);
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - delQs;
            elseif fv.splines(kVV,2).v2 < Vm && Vm <= fv.splines(kVV,3).v1
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            elseif fv.splines(kVV,3).v1 < Vm && Vm <= fv.splines(kVV,3).v2
                delQs = del(fv.splines(kVV,3),Vm);
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV) = dQkm_dVm(i) - delQs;
            elseif fv.splines(kVV,3).v2 < Vm && Vm <= fv.splines(kVV,4).v1
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - ((Q4-Q3)/(V4-V3));
            elseif fv.splines(kVV,4).v1 < Vm && Vm <= fv.splines(kVV,4).v2
                delQs = del(fv.splines(kVV,4),Vm);
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - delQs;
            elseif fv.splines(kVV,4).v2 < Vm 
                Y(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            end
            
            kVV = kVV + 1;
        end
    end

    
    % Z (4Nfv+Nfv+2Nfv,fpf+4Nfv,vw x 4Nfv+Nfv+2Nfv,fpf+4Nfv,vw) 
    %           delG / delVarInternal
    Z = zeros(5*Nfv+2*fv.Nfpf+4*fv.Nvw, 5*Nfv+2*fv.Nfpf+4*fv.Nvw);
    
    kDeltaPhi1 = 1;         
    kFPF=1; kVW=1; kVV=1;
    for i=1:Nfv

        % deltaG1/Vfv
        Z(i, i) = - fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i)))/(fv.a(i)*fv.Vt(i));
        % deltaG1/Ifv
        Z(i, Nfv + i) = - 1;

        % deltaG3/Vfv
        Z(Nfv + i, i) = fv.eta(i)*fv.Ifv(i) - dPkm_dVfv(i);
        % deltaG3/Ifv
        Z(Nfv + i, Nfv + i) = fv.eta(i)*fv.Vfv(i);
        % deltaG3/m
        Z(Nfv + i, 2*Nfv + i) = - dPkm_dm(i);
        % deltaG3/Alpha
        Z(Nfv + i, 3*Nfv + i) = - dPkm_dalpha(i);
        
        % deltaPhi2/Vfv
        Z(2*Nfv+fv.Nfpf+fv.Nvv + i, i) = FunDelFBrelax((-fv.eta(i)*fv.Ifv(i)), (fv.Snom(i)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.xiC(i), xi, 1);
%         ( -fv.eta(i)*fv.Ifv(i) ) * ...
%             ( (fv.Snom(i) - fv.eta(i)*fv.Vfv(i)*fv.Ifv(i))/(sqrt((fv.Snom(i) - fv.eta(i)*fv.Vfv(i)*fv.Ifv(i))^2 + fv.xiC(i)^2 + xi)) - 1);
        % deltaPhi2/Ifv
        Z(2*Nfv+fv.Nfpf+fv.Nvv + i, Nfv + i) = FunDelFBrelax((-fv.eta(i)*fv.Vfv(i)), (fv.Snom(i)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.xiC(i), xi, 1);
%         ( -fv.eta(i)*fv.Vfv(i) ) * ...
%             ( (fv.Snom(i) - fv.eta(i)*fv.Vfv(i)*fv.Ifv(i))/(sqrt((fv.Snom(i) - fv.eta(i)*fv.Vfv(i)*fv.Ifv(i))^2 + fv.xiC(i)^2 + xi)) - 1);
        % deltaPhi2/xiC
        Z(2*Nfv+fv.Nfpf+fv.Nvv + i, 4*Nfv + i) = FunDelFBrelax(0, (fv.Snom(i)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.xiC(i), xi, 2);
%         ( fv.xiC(i)/(sqrt((fv.Snom(i) - fv.eta(i)*fv.Vfv(i)*fv.Ifv(i))^2 + fv.xiC(i)^2 + xi)) ) - 1;

        if fv.cont(i) == 1              %  FPF
            k_aux = tan(acos(fv.FPF(2)));
            
            % deltaPhi1/Vfv
            Z(2*Nfv + kDeltaPhi1, i) = - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) / ((fv.a(i)*fv.Vt(i))^2);
            % deltaPhi1/Ifv
            Z(2*Nfv + kDeltaPhi1, Nfv + i) = 1;
            % deltaPhi1/xiC
            Z(2*Nfv + kDeltaPhi1, 4*Nfv + i) = 1;
            
            % deltaPhiFPF1/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, i) = dQkm_dVfv(i) - dPkm_dVfv(i)*k_aux;
            % deltaPhiFPF1/m
            Z(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, 2*Nfv + i) = dQkm_dm(i) - dPkm_dm(i)*k_aux;
            % deltaPhiFPF1/Alpha
            Z(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, 3*Nfv + i) = dQkm_dalpha(i) - dPkm_dalpha(i)*k_aux;
            % deltaPhiFPF1/Qli
            Z(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, 5*Nfv + kFPF) = -1;
            % deltaPhiFPF1/Qls
            Z(3*Nfv+fv.Nfpf+fv.Nvv + kFPF, 5*Nfv+fv.Nfpf + kFPF) = 1;
            
            % deltaPhiFPF2/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+fv.Nfpf + kFPF, i) = FunDelFBrelax((dQkm_dVfv(i)+dQdisp_dVfv(i)), (fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi, 1); 
            % deltaPhiFPF2/m
            Z(3*Nfv+fv.Nfpf+fv.Nvv+fv.Nfpf + kFPF, 2*Nfv + i) = FunDelFBrelax((dQkm_dm(i)+dQdisp_dm(i)), (fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi, 1); 
            % deltaPhiFPF2/Alpha
            Z(3*Nfv+fv.Nfpf+fv.Nvv+fv.Nfpf + kFPF, 3*Nfv + i) = FunDelFBrelax((dQkm_dalpha(i)+dQdisp_dalpha(i)), (fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi, 1); 
            % deltaPhiFPF2/Qli
            Z(3*Nfv+fv.Nfpf+fv.Nvv+fv.Nfpf + kFPF, 5*Nfv + kFPF) = FunDelFBrelax(0, (fv.Qkm(i)+fv.Qdisp(i)), fv.Qli(kFPF), xi, 2); 
            
            % deltaPhiFPF3/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+2*fv.Nfpf + kFPF, i) = FunDelFBrelax((dQdisp_dVfv(i)-dQkm_dVfv(i)), (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi, 1);
            % deltaPhiFPF3/m
            Z(3*Nfv+fv.Nfpf+fv.Nvv+2*fv.Nfpf + kFPF, 2*Nfv + i) = FunDelFBrelax((dQdisp_dm(i)-dQkm_dm(i)), (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi, 1);
            % deltaPhiFPF3/Alpha
            Z(3*Nfv+fv.Nfpf+fv.Nvv+2*fv.Nfpf + kFPF, 3*Nfv + i) = FunDelFBrelax((dQdisp_dalpha(i)-dQkm_dalpha(i)), (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi, 1);
            % deltaPhiFPF3/Qls
            Z(3*Nfv+fv.Nfpf+fv.Nvv+2*fv.Nfpf + kFPF, 5*Nfv+fv.Nfpf + kFPF) = FunDelFBrelax(0, (fv.Qdisp(i)-fv.Qkm(i)), fv.Qls(kFPF), xi, 2);
            
            kFPF = kFPF + 1;
            kDeltaPhi1 = kDeltaPhi1 + 1;
            
        elseif fv.cont(i) == 2         %  VW     
            P1 = fv.VW(3)*fv.Snom(i); P2 = fv.VW(5)*fv.Snom(i);
            
            % deltaPhiVW1/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf + kVW, i) = - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) / ((fv.a(i)*fv.Vt(i))^2);
            % deltaPhiVW1/Ifv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf + kVW, Nfv + i) = 1;
            % deltaPhiVW1/xiC
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf + kVW, 4*Nfv + i) = 1;
            % deltaPhiVW1/Pls
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf + kVW, 5*Nfv+2*fv.Nfpf+fv.Nvw + kVW) = 1;
            
            % deltaPhiVW2/Pmax
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf + kVW) = 1;
            % deltaPhiVW2/PmaxLI
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf+2*fv.Nvw + kVW) = -1;
            % deltaPhiVW2/PmaxLS
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf+3*fv.Nvw + kVW)= 1;
            
            % deltaPhiVW3/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+2*fv.Nvw + kVW, i) = FunDelFBrelax((-fv.eta(i)*fv.Ifv(i)), (fv.Pmax(kVW)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.Pls(kVW), xi, 1);
            % deltaPhiVW3/Ifv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+2*fv.Nvw + kVW, Nfv + i) = FunDelFBrelax((-fv.eta(i)*fv.Vfv(i)), (fv.Pmax(kVW)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.Pls(kVW), xi, 1);
            % deltaPhiVW3/Pmax
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+2*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf + kVW) = FunDelFBrelax(1, (fv.Pmax(kVW)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.Pls(kVW), xi, 1);
            % deltaPhiVW3/Pls
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+2*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf+fv.Nvw + kVW) = FunDelFBrelax(0, (fv.Pmax(kVW)-fv.eta(i)*fv.Vfv(i)*fv.Ifv(i)), fv.Pls(kVW), xi, 2);
            
            % deltaPhiVW4/Pmax
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+3*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf + kVW) = FunDelFBrelax(-1, (P1-fv.Pmax(kVW)), fv.PmaxLS(kVW), xi, 1);
            % deltaPhiVW4/PmaxlS
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+3*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf+3*fv.Nvw + kVW) = FunDelFBrelax(0, (P1-fv.Pmax(kVW)), fv.PmaxLS(kVW), xi, 2);
            
            % deltaPhiVW5/Pmax
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+4*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf + kVW) = FunDelFBrelax(1, (fv.Pmax(kVW)-P2), fv.PmaxLI(kVW), xi, 1);
            % deltaPhiVW5/PmaxLI
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+4*fv.Nvw + kVW, 5*Nfv+2*fv.Nfpf+2*fv.Nvw + kVW) = FunDelFBrelax(0, (fv.Pmax(kVW)-P2), fv.PmaxLI(kVW), xi, 2);
            
            % deltaPhiVW6/Vfv
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+5*fv.Nvw + kVW, i) = dQkm_dVfv(i);
            % deltaPhiVW6/m
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+5*fv.Nvw + kVW, 2*Nfv + i) = dQkm_dm(i);
            % deltaPhiVW6/Alpha
            Z(3*Nfv+fv.Nfpf+fv.Nvv+3*fv.Nfpf+5*fv.Nvw + kVW, 3*Nfv + i) = dQkm_dalpha(i);
            
            kVW = kVW + 1;
            
        elseif fv.cont(i) == 3          %  VV
            
            % deltaPhi1/Vfv
            Z(2*Nfv + kDeltaPhi1, i) = - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) / ((fv.a(i)*fv.Vt(i))^2);
            % deltaPhi1/Ifv
            Z(2*Nfv + kDeltaPhi1, Nfv + i) = 1;
            % deltaPhi1/xiC
            Z(2*Nfv + kDeltaPhi1, 4*Nfv + i) = 1;
            
            % deltaPhiVV/Vfv 
            Z(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, i) = dQkm_dVfv(i);
            % deltaPhiVV/m
            Z(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, 2*Nfv + i) = dQkm_dm(i);
            % deltaPhiVV/Alpha
            Z(3*Nfv+4*fv.Nfpf+fv.Nvv+6*fv.Nvw + kVV, 3*Nfv + i) = dQkm_dalpha(i);
            
            kVV = kVV + 1;
            kDeltaPhi1 = kDeltaPhi1 + 1;
            
        end
    end
end


J = [Jgrid X; Y Z];
    
end