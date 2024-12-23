%================ POWER CALCULATION ======================
% Of every buses, including slack

function [Pcalc,Qcalc,soma1,soma3] = calcPot(DBAR,DRAM,Nb,Nramo)

a = ones(Nramo,1)./DRAM.a;      
Pcalc = zeros(Nb,1);
Qcalc = zeros(Nb,1);
soma1 = zeros(Nb,1);  soma3 = zeros(Nb,1);      

for k=1:Nb             
    soma2=0; soma4=0;
    for m=1:Nb                                  % Search for neighbours
        if m==k                                 
            continue 
        end        
        [km,flag] = km_flag(Nramo, DRAM, k, m);     
        if km~=0                               
            if flag==1                        
                soma1(k) = soma1(k) + a(km)^2*DRAM.g(km);
                soma2 = soma2 + DBAR.V(m)*a(km)*(DRAM.g(km)*cos(DBAR.theta(k) - DBAR.theta(m)) + ...
                    DRAM.b(km)*sin(DBAR.theta(k) - DBAR.theta(m)));
                soma3(k) = soma3(k) + DRAM.Bsh(km) + a(km)^2*DRAM.b(km);
                soma4 = soma4 + DBAR.V(m)*a(km)*(DRAM.g(km)*sin(DBAR.theta(k) - DBAR.theta(m)) - ...
                    DRAM.b(km)*cos(DBAR.theta(k) - DBAR.theta(m)));
            elseif flag==0                     
                soma1(k) = soma1(k) + DRAM.g(km);
                soma2 = soma2 + DBAR.V(m)*a(km)*(DRAM.g(km)*cos(DBAR.theta(k) - DBAR.theta(m)) + ...
                    DRAM.b(km)*sin(DBAR.theta(k) - DBAR.theta(m)));
                soma3(k) = soma3(k) + DRAM.Bsh(km) + DRAM.b(km);
                soma4 = soma4 + DBAR.V(m)*a(km)*(DRAM.g(km)*sin(DBAR.theta(k) - DBAR.theta(m)) - ...
                    DRAM.b(km)*cos(DBAR.theta(k) - DBAR.theta(m)));
            end
        end
    end
    Pcalc(k) = DBAR.V(k)^2*soma1(k) - DBAR.V(k)*soma2;
    Qcalc(k) = -(DBAR.V(k)^2*(DBAR.mhoSh(k) + soma3(k))) - DBAR.V(k)*soma4;  
end
