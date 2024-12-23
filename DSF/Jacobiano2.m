%% ========== JACOBIAN =================

function [J] = Jacobiano2(DBAR,DRAM,Nb,Nramo,Pcalc,Qcalc,som1,som3)

a = ones(Nramo,1)./DRAM.a;      

H = zeros(DBAR.NPQ+DBAR.NPV, DBAR.NPQ+DBAR.NPV);        % dDeltaP/dtheta
N = zeros(DBAR.NPQ+DBAR.NPV, DBAR.NPQ);                 % dDeltaP/DV
M = zeros(DBAR.NPQ, DBAR.NPQ+DBAR.NPV);                 % dDeltaQ/dtheta
L = zeros(DBAR.NPQ, DBAR.NPQ);                          % dDeltaQ/dV

% == H - delDeltaP/delTheta (NPV+NPQ x NPV+NPQ) (Nb-1 x Nb-1)
for i=2:Nb
    for j=2:Nb
        [km,~] = km_flag(Nramo, DRAM, i, j);
        if i==j        
            H(i-1, j-1) = H(i-1, j-1) + Qcalc(i) + (DBAR.V(i)^2* (abs(DBAR.mhoSh(i)) + som3(i))); 
        else            
            if km~=0    
                H(i-1, j-1) = + DBAR.V(i)*DBAR.V(j)*a(km) * (DRAM.g(km)*sin(DBAR.theta(i) - DBAR.theta(j)) - ...
                    DRAM.b(km)*cos(DBAR.theta(i) - DBAR.theta(j)));                  
            end
        end
    end
end

% == N - delDeltaP/delV (NPV+NPQ x NPQ)
for i=2:Nb
    for j=1:DBAR.NPQ       % DBAR.NFV is the number of bus with PV system
        n = DBAR.index_PQ(j);
        [km,~] = km_flag(Nramo, DRAM, i, n);
        
        if i==n   
            N(i-1, j) = N(i-1 ,j) - (Pcalc(i) + DBAR.V(i)^2*som1(i))/DBAR.V(i);
        else        
            if km~=0                
                N(i-1, j) = + DBAR.V(i)*a(km) * (DRAM.g(km)*cos(DBAR.theta(i) - DBAR.theta(n)) + ...
                    DRAM.b(km)*sin(DBAR.theta(i) - DBAR.theta(n)));
            end
        end
    end
end


% == M  - delDeltaQ/delTheta (NPQ x NPV+NPQ)
for i=1:DBAR.NPQ
    for j=2:Nb
        n = DBAR.index_PQ(i);
        [km,~] = km_flag(Nramo, DRAM, n, j);
        if n==j    
            M(i, j-1) = M(i, j-1) - Pcalc(n) + (DBAR.V(n)^2*som1(n));
        else       
            if km~=0
                M(i,j-1) = - DBAR.V(n)*DBAR.V(j)*a(km)*(DRAM.g(km)*cos(DBAR.theta(n) - DBAR.theta(j)) + ...
                    DRAM.b(km)*sin(DBAR.theta(n) - DBAR.theta(j)));
            end
        end
    end
end

% == L  - delDeltaQ/delV (NPQ x NPQ)
for i=1:DBAR.NPQ
    n = DBAR.index_PQ(i);
    for j=1:DBAR.NPQ
        k = DBAR.index_PQ(j);
        [km,~] = km_flag(Nramo, DRAM, n, k);
        if i==j          
            L(i,j) = L(i,j) - (Qcalc(n) - DBAR.V(n)^2*(abs(DBAR.mhoSh(n)) + som3(n))) / DBAR.V(n);
        else                
            if km~=0
                L(i,j) = + DBAR.V(n)*a(km) * (DRAM.g(km)*sin(DBAR.theta(n) - DBAR.theta(k)) - ...
                    DRAM.b(km)*cos(DBAR.theta(n) - DBAR.theta(k)));
            end            
        end
    end
end

J = [H N; M L];

end
