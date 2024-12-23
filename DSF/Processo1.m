%% ====== SCPRIT PROCESS 1 =========
% Calculation of PVs power injections 
% Update Pesp and Qesp of PV system buses

for i=1:Nfv
    
    % -- Pcc
    PTvalue = PTcurve(fv.PTcurve(i,2:end),fv.T(i));         % func calc PTvalue via piece-wise linear curve
    fv.Pcc(i) = fv.Pnom(i)*fv.G(i)*PTvalue;
    
    % -- Active power AC 
    % limit
    fv.Plimit_min(i) = min([fv.Snom(i) , fv.Plimit(i)]);  % Plimit defined by VW
    
    if fv.Pcc(i)*fv.eta(i) >= fv.Plimit_min(i)
        fv.Pca_lin(i) = fv.Plimit_min(i);
    else
        fv.Pca_lin(i) = fv.Pcc(i)*fv.eta(i);
    end
    
    % -- Reactive power 
    if fv.cont(i) == -1                             % FPF 
        fv.Qca_lin(i) = sign(fv.fp(i))*fv.Pca_lin(i)*tan(acos(abs(fv.fp(i))));
    elseif fv.cont(i) == 1 || fv.cont(i) == 3       % VV 
        fv.Qca_lin(i) = fv.Qdesired(i);     % defined in Process 2 (invcontrol)
    end
    
    % -- Verify inverter capacity
    if (fv.Pca_lin(i)^2 + fv.Qca_lin(i)^2) <= fv.Snom(i)^2      % OK
        fv.Pca(i) = fv.Pca_lin(i);     
        fv.Qca(i) = fv.Qca_lin(i); 
    elseif (fv.Pca_lin(i)^2 + fv.Qca_lin(i)^2) > fv.Snom(i)^2   % violation
        % priorities
        if fv.fpPri(i) == 1        % priority FP
            fv.Pca(i) = fv.Snom(i)*abs(fv.fp(i));
            fv.Qca(i) = fv.Snom(i)*sqrt(1 - fv.fp(i)^2)*sign(fv.fp(i));
        elseif fv.fpPri(i) == 0 && fv.wattPri(i) == 1       % prioritt watt
            if fv.Pca_lin(i) >= fv.Snom(i)
                fv.Pca(i) = fv.Snom(i);
            else 
                fv.Pca(i) = fv.Pca_lin(i);
            end
            fv.Qca(i) = sqrt(fv.Snom(i)^2 - fv.Pca(i)^2) * sign(fv.Qca_lin(i));
        elseif fv.fpPri(i) == 0 && fv.wattPri(i) == 0       % priority var
            if fv.Qca_lin(i) >= fv.Snom(i)
                fv.Qca(i) = fv.Snom(i);
            else
                fv.Qca(i) = fv.Qca_lin(i);
            end
            fv.Pca(i) = sqrt(fv.Snom(i)^2 - fv.Qca(i)^2);
        end
    end
        
    
    % --- Update injection power 
    barFV = fv.bar(i);
    DBAR.Pgen(barFV) = fv.Pca(i);
    DBAR.Qgen(barFV) = fv.Qca(i);
    DBAR.Pesp(barFV) = DBAR.Pgen(barFV) - DBAR.Pc(barFV);
    DBAR.Qesp(barFV) = DBAR.Qgen(barFV) - DBAR.Qc(barFV);
    
    
end