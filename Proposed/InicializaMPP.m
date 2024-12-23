function [Vmpp,Impp] = InicializaMPP(fv,Nfv)
% find the estimate MPP for every PV system
% OUTPUT - Vmpp (Nfv x 1)
%        - Impp (Nfv x 1)

X = zeros(101,1); I = zeros(101,1);
Vmpp = zeros(Nfv,1); Impp = zeros(Nfv,1);

for i=1:Nfv
    k=1;
    for V=0:(fv.Voc(i)*fv.Nss(i)/fv.Vbase(i))/100:(fv.Voc(i)*fv.Nss(i)/fv.Vbase(i))
        
        I(k) = fv.Iph(i) - fv.Is(i)*(exp(V/(fv.a(i)*fv.Vt(i))) - 1);
        X(k) = V;
        k = k + 1;
    end
    P = X.*I;
    [~,endMPP] = max(P);
    Vmpp(i) = X(endMPP);
    Impp(i) = I(endMPP);    
end
