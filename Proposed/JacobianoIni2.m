%% ================== JACOBIAN ===========================

function [J] = JacobianoIni2(DBAR,DRAM,fv,Nb,Nramo,Nfv,Pcalc,Qcalc,som1,som3)
%% === Constant

a = ones(Nramo,1)./DRAM.a;    


%% == deltaSkm
if ~isempty(fv)
    
    dPkm_dTm = zeros(Nfv,1);  dPkm_dVm = zeros(Nfv,1);
    dPkm_dVfv = zeros(Nfv,1); 
    dPkm_dm = zeros(Nfv,1); 
    dPkm_dalpha = zeros(Nfv,1);
    
    dQkm_dTm = zeros(Nfv,1);  dQkm_dVm = zeros(Nfv,1);
    dQkm_dVfv = zeros(Nfv,1);
    dQkm_dm = zeros(Nfv,1); dQkm_dalpha = zeros(Nfv,1);
    
%     dPmk_dTm = zeros(Nfv,1); dPmk_dVm = zeros(Nfv,1);
%     dPmk_dVfv = zeros(Nfv,1); dPmk_dm = zeros(Nfv,1);
%     dPmk_dalpha = zeros(Nfv,1);
%     dQmk_dTm = zeros(Nfv,1); dQmk_dVm = zeros(Nfv,1);
%     dQmk_dVfv = zeros(Nfv,1); dQmk_dm = zeros(Nfv,1);
%     dQmk_dalpha = zeros(Nfv,1);
    
%     dQdisp_dTm = zeros(Nfv,1); dQdisp_dVm = zeros(Nfv,1);
%     dQdisp_dVfv = zeros(Nfv,1);
%     dQdisp_dm = zeros(Nfv,1); dQdisp_dalpha = zeros(Nfv,1);

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

%         % -- dPmk
%         dPmk_dTm(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
%             ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
%         dPmk_dVm(i) = 2*DBAR.V(barFV)*fv.g_t(i) - sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
%             ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
%         dPmk_dVfv(i) = - sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...
%             ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
%         dPmk_dm(i) = - sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...
%             ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
%         dPmk_dalpha(i) = - sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
%             ( fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) - fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
%         
%         % -- dQmk
%         dQmk_dTm(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
%             ( - fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) - fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
%         dQmk_dVm(i) = -2*DBAR.V(barFV)*fv.b_t(i) + sqrt(3/8)*fv.m(i)*fv.Vfv(i) * ...
%             ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
%         dQmk_dVfv(i) = sqrt(3/8)*fv.m(i)*DBAR.V(barFV) * ...                       
%             ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
%         dQmk_dm(i) = sqrt(3/8)*fv.Vfv(i)*DBAR.V(barFV) * ...                       
%             ( - fv.g_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) );
%         dQmk_dalpha(i) = sqrt(3/8)*fv.m(i)*fv.Vfv(i)*DBAR.V(barFV) * ...
%             ( fv.g_t(i)*cos(DBAR.theta(barFV) - fv.alpha(i)) + fv.b_t(i)*sin(DBAR.theta(barFV) - fv.alpha(i)) );
        
%         % -- Qdisp
%         dQdisp_dTm(i) = - dPkm_dTm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
%         dQdisp_dVm(i) = - dPkm_dVm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
%         dQdisp_dVfv(i) = - dPkm_dVfv(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
%         dQdisp_dm(i) = - dPkm_dm(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
%         dQdisp_dalpha(i) = - dPkm_dalpha(i) * ( fv.Pkm(i)/fv.Qdisp(i) );
            
    end
end

%%  Grid

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
                    H(i-1, j-1) = H(i-1, j-1) + dPkm_dTm(m);  % adiciona a derivada do Pkm
%                     H(i-1, j-1) = H(i-1, j-1) - dPmk_dTm(m);  % subtrai a derivada do Pmk
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


%% ==== jacobian extension

X = [];
Y = [];
Z = [];

if ~isempty(fv)    
    
    % X (2NPQ+2NFV+NPV  X  4Nfv) 
    %       delDeltaP e delDeltaQ / delVarInternals
    X = zeros(2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV, 4*Nfv);

    for i=1:Nfv
        barFV = fv.bar(i);
        indexFV = find(auxIndex_PQ_FV==barFV);  % bus PV system postiion in deltaQ
        
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

    
    % Y (4Nfv X 2NPQ+2Nfv+NPV) 
    %       delG/delTheta e delG/delVgrid
    Y = zeros(4*Nfv, 2*(DBAR.NPQ+DBAR.NFV)+DBAR.NPV);
    
    kFPF=1; kVW=1; kVV=1;
    for i=1:Nfv
        barFV = fv.bar(i);
        indexFV = find(auxIndex_PQ_FV==barFV); 

        % deltaG3/Theta_m
        Y(2*Nfv + i, barFV-1) = -dPkm_dTm(i);
        % deltaG3/Vm
        Y(2*Nfv + i, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = -dPkm_dVm(i);
        
        
        if fv.cont(i) == 1          % control FPF      
            % deltaGfpf/Theta_m
            Y(3*Nfv + kFPF, barFV-1) = dQkm_dTm(i) - dPkm_dTm(i)*tan(acos(fv.FPF(2)));
            % deltaGfpf/Vm
            Y(3*Nfv + kFPF, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i) - dPkm_dVm(i)*tan(acos(fv.FPF(2)));
            
            kFPF = kFPF+1;
            
        elseif fv.cont(i) == 2      % control VW            
            % deltaGvw/Theta_m
            Y(3*Nfv+fv.Nfpf + kVW, barFV-1) = dQkm_dTm(i);
            % deltaGvw/Vm
            Y(3*Nfv+fv.Nfpf + kVW, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            
            kVW = kVW + 1;
            
        elseif fv.cont(i) == 3      % control VV
            
            % deltaGvv/Theta_m
            Y(3*Nfv+fv.Nfpf+fv.Nvw + kVV, barFV-1) = dQkm_dTm(i);
            % deltaGvv/Vm
            Y(3*Nfv+fv.Nfpf+fv.Nvw + kVV, DBAR.NPQ+DBAR.NFV+DBAR.NPV + indexFV) = dQkm_dVm(i);
            
            kVV = kVV + 1;
        end
    end

    
    % Z (4Nfv x 4Nfv) 
    %           delG / delVarInternals
    Z = zeros(4*Nfv, 4*Nfv);
    
    kFPF=1; kVW=1; kVV=1;
    for i=1:Nfv

        % deltaG1/Vfv
        Z(i, i) = - fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i)))/(fv.a(i)*fv.Vt(i));
        % deltaG1/Ifv
        Z(i, Nfv + i) = - 1;
        
        % deltaG2/Vfv
        Z(Nfv + i, i) = - fv.Vfv(i)*fv.Is(i)*exp(fv.Vfv(i)/(fv.a(i)*fv.Vt(i))) / ((fv.a(i)*fv.Vt(i))^2);
        % deltaG2/Ifv
        Z(Nfv + i, Nfv + i) = 1;
        
        % deltaG3/Vfv
        Z(2*Nfv + i, i) = fv.eta(i)*fv.Ifv(i) - dPkm_dVfv(i);
        % deltaG3/Ifv
        Z(2*Nfv + i, Nfv + i) = fv.eta(i)*fv.Vfv(i);
        % deltaG3/m
        Z(2*Nfv + i, 2*Nfv + i) = - dPkm_dm(i);
        % deltaG3/Alpha
        Z(2*Nfv + i, 3*Nfv + i) = - dPkm_dalpha(i);

        if fv.cont(i) == 1              % control FPF
            k_aux = tan(acos(fv.FPF(2)));
            
            % deltaGfpf/Vfv
            Z(3*Nfv + kFPF, i) = dQkm_dVfv(i) - dPkm_dVfv(i)*k_aux;
            % deltaGfpf/m
            Z(3*Nfv + kFPF, 2*Nfv + i) = dQkm_dm(i) - dPkm_dm(i)*k_aux;
            % deltaGfpf/Alpha
            Z(3*Nfv + kFPF, 3*Nfv + i) = dQkm_dalpha(i) - dPkm_dalpha(i)*k_aux;
            
            kFPF = kFPF + 1;
            
        elseif fv.cont(i) == 2         % control VW                 
            % deltaGvw/Vfv
            Z(3*Nfv+fv.Nfpf + kVW, i) = dQkm_dVfv(i);
            % deltaGvw/m
            Z(3*Nfv+fv.Nfpf + kVW, 2*Nfv + i) = dQkm_dm(i);
            % deltaGvw/Alpha
            Z(3*Nfv+fv.Nfpf + kVW, 3*Nfv + i) = dQkm_dalpha(i);
            
            kVW = kVW + 1;
            
        elseif fv.cont(i) == 3          % control VV
            % deltaPhiVV/Vfv 
            Z(3*Nfv+fv.Nfpf+fv.Nvw + kVV, i) = dQkm_dVfv(i);
            % deltaPhiVV/m
            Z(3*Nfv+fv.Nfpf+fv.Nvw + kVV, 2*Nfv + i) = dQkm_dm(i);
            % deltaPhiVV/Alpha
            Z(3*Nfv+fv.Nfpf+fv.Nvw + kVV, 3*Nfv + i) = dQkm_dalpha(i);
            
            kVV = kVV + 1;
            
        end
    end
end


J = [Jgrid X; Y Z];
    
end