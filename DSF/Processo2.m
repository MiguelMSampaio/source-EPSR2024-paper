%% ========== SCRIPT PROCESS 2 - CONTROL ACTIONS =======
% for the PV systems with invcontrol, update Qdesired e Plimit

for i = 1:fv.NinvCtrl
    fvI = fv.index_invCtrl(i);         
    barI = fv.bar(fvI);              
    
    refVar = fv.invControl(fvI,1);    
    deltaQ = fv.invControl(fvI,2);
    
    refW = fv.invControl(fvI,11);
    deltaP = fv.invControl(fvI,12);
    
    % ====== Calculation of base power
    % Active
    if refW == 1        % Pnom
        fv.Pbase(i) = fv.Pnom(fvI);
    elseif refW == 2    % available
        fv.Pbase(i) = fv.Pcc(fvI)*fv.eta(fvI);
    elseif refW == 3    % Snom
        fv.Pbase(i) = fv.Snom(fvI);
    end
    % Reactive
    if refVar == 1      % Qdisp
        fv.Qbase(i) = sqrt(fv.Snom(fvI)^2 - fv.Pca(fvI)^2);
    elseif refVar == 2  % Snom
        fv.Qbase(i) = fv.Snom(fvI);
    end
    
    % ===== VOLT-VAR
    if fv.cont(fvI) == 1
        
        % - look-up table VV
        fv.Qd_fun(i) = FunVV(fv.Vmon(i,itCtrl+1), fv.invControl(fvI,3:10), fv.Qbase(i));
        
        % - Define inverter capacity limits
        if refVar == 2 && fv.wattPri(fvI) == 1      % axis Y pu of Snom and priority watt
            fv.Qd_lim(i) = min( sqrt(fv.Snom(fvI)^2 - fv.Pca(fvI)^2), fv.Snom(fvI) );
        else
            fv.Qd_lim(i) = fv.Snom(fvI);
        end
        
        % - Final values
        fv.Qd_end(i) = min( abs(fv.Qd_fun(i)), abs(fv.Qd_lim(i)) ) * sign(fv.Qd_fun(i));
        
        % - Update scale factor
        if fv.FdeltaQ(i) == -1
            fv.invControl(fvI,2) = AtualizaFE(deltaQ,passoV(i),passoVold(i));
            deltaQ = fv.invControl(fvI,2);
        end
                
        % - Update Qdesired
        fv.Qdesired(fvI) = fv.Qca(fvI) + (fv.Qd_end(i) - fv.Qca(fvI)) * deltaQ;
        
    % ===== VOLT-WATT
    elseif fv.cont(fvI) == 2 
        
        % - Look-up table
        fv.Pl_fun(i) = FunVW(fv.Vmon(i,itCtrl+1), fv.invControl(fvI,13:16), fv.Pbase(i));
        
%         % - Define inverter capacity limits
%         fv.Pl_lim(i) = min( sqrt(fv.Snom(fvI)^2 - fv.Qca(fvI)^2), fv.Pnom(fvI) );     
%         
%         % - Final values
%         fv.Pl_end(i) = min(fv.Pl_fun(i), fv.Pl_lim(i));
        fv.Pl_end(i) = fv.Pl_fun(i);
        
        % - Update scale factor
        if fv.FdeltaP(i) == -1
            fv.invControl(fvI,12) = AtualizaFE(deltaP,passoV(i),passoVold(i));
            passoP = fv.invControl(fvI,12);
        end
        
        % - Update Qdesired
        fv.Plimit(fvI) = fv.Pca(fvI) + (fv.Pl_end(i) - fv.Pca(fvI))*deltaP;
        
    % ===== VV+VW
    elseif fv.cont(fvI) == 3
        
         % - look-up table VV
        fv.Qd_fun(i) = FunVV(fv.Vmon(i,itCtrl+1), fv.invControl(fvI,3:10), fv.Qbase(i));
        fv.Pl_fun(i) = FunVW(fv.Vmon(i,itCtrl+1), fv.invControl(fvI,15:18), fv.Pbase(i));
        
        % - Define inverter capacity limit
        if refVar == 2 && fv.wattPri(fvI) == 1      % axis Y pu of Snom and priority Watt 
            fv.Qd_lim(i) = min( sqrt(fv.Snom(fvI)^2 - fv.Pca(fvI)^2), fv.Snom(fvI) );
        else
            fv.Qd_lim(i) = fv.Snom(fvI);
        end
        fv.Pl_fun(i) = FunVW(fv.Vmon(i,itCtrl+1), fv.invControl(fvI,15:18), fv.Pbase(i));
        
        % - Final values
        fv.Qd_end(i) = min( abs(fv.Qd_fun(i)), abs(fv.Qd_lim(i)) ) * sign(fv.Qd_fun(i));
        fv.Pl_end(i) = min(fv.Pl_fun(i), fv.Pl_lim(i));
        
        % - Update scale factor
        if fv.FdeltaQ(i) == -1
            fv.invControl(fvI,2) = AtualizaFE(deltaQ,passoV(i),passoVold(i));
            deltaQ = fv.invControl(fvI,2);
            fv.invControl(fvI,12) = AtualizaFE(deltaP,passoV(i),passoVold(i));
            passoP = fv.invControl(fvI,12);
        end
                
        % - Update Qdesired
        fv.Qdesired(fvI) = fv.Qca(fvI) + (fv.Qd_end(i) - fv.Qca(fvI)) * deltaQ;
        fv.Plimit(fvI) = fv.Pca(fvI) + (fv.Pl_end(i) - fv.Pca(fvI))*deltaP;
        
    end
    
end

