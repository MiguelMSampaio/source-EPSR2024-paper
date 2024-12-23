%% ====== FUNC TO UPDATE SCALE FACTOR ======
% based on monitored voltage step and last it ctrl 
% update every scale factor individualy

function FE = AtualizaFE(delta,passoV,passoVold)
%   ENTRADA         - delta: actual scale factor
%                   - passoV: actual monitored voltage step
%                   - passoVold: last it monitored voltage step

%   SAIDA           - FE: new scale factor

FE = delta;
% fprintf('%.3f\n',passoV/passoVold)

if passoVold ~= 0
    if abs(passoV) > 0.8*abs(passoVold) && delta > 0.2
        FE = delta - 0.1;
    elseif abs(passoV) > 0.6*abs(passoVold) && delta > 0.2
        FE = delta - 0.05;
    elseif abs(passoV) < 0.2*abs(passoVold) && delta < 0.9
        FE = delta + 0.1;
    elseif abs(passoV) < 0.4*abs(passoVold) && delta < 0.9
        FE = delta + 0.05;
    end
end

end