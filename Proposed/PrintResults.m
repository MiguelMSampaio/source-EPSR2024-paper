%% ==== FUNC TO PRINT RESULTS ===== 

function PrintResults(ps)

fprintf('  =================  Fluxo de Potência %s - Método Newton Raphson  =================== \n',ps.sistema);
fprintf('  Nº de iterações total     : %d \n', ps.results.gerais(1));
fprintf('  Nº de iterações parcial   : %d \n', ps.results.gerais(7));
fprintf('  Resíduo máximo            : %g \n', ps.results.gerais(2));
% fprintf('  Tempo de execução total   : %.4f s\n', ps.results.gerais(3));
fprintf('  Tempo de execução método  : %.4f s\n', ps.results.gerais(4));

% fprintf('\n----------- Variáveis de estado\n'); 
% fprintf('%5s %10s %10s %10s %10s %10s %10s\n','Barra','Tensão(pu)','Theta(deg)','Pgen(MW)','Qgen(MVAr)','Pc(MW)','Qc(MVAr)')
% fprintf('%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.barra.')

% if ps.results.fv.Nfv>0
% %     fprintf('\n----------- Variáveis de internas\n'); 
% %     fprintf('%5s %10s %10s %10s %10s %13s %10s %10s %10s %10s %10s\n','BarFV','Vfv(pu)','Ifv(pu)','Pfv(pu)','m',...
% %         'alpha(deg)','Vk(pu)','Pkm(pu)','Qkm(pu)','Skm(pu)','fp');
% %     fprintf('%5d %10.4f %10.4f %10.4f %10.4f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.fv.varInt.');
%     
%     ps.results.fv.varInt(:,[4 8:10]) = ps.results.fv.varInt(:,[4 8:10])*ps.Sbase/1e6;
%     fprintf('\n----------- Variáveis de internas\n'); 
%     fprintf('%5s %10s %10s %10s %10s %13s %10s %10s %10s %10s %10s\n','BarFV','Vfv(pu)','Ifv(pu)','Pfv(MW)','m',...
%         'alpha(deg)','Vk(pu)','Pkm(MW)','Qkm(MW)','Skm(MW)','fp');
%     fprintf('%5d %10.4f %10.4f %10.4f %10.4f %13.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.fv.varInt.');
%     
%     fprintf('\n        --- Variáveis de fogla FPF\n')
%     fprintf('%5s %10s %10s %10s\n', 'BarFV','xiC','Qli','Qls')
%     fprintf('%5d %10.5f %10.5f %10.5f\n',ps.results.fv.varFolgaFPF.')
%     
%     fprintf('\n        --- Variáveis de fogla VW\n')
%     fprintf('%5s %10s %10s %10s %10s %10s \n', 'BarFV','xiC','Pmax','Pls','PmaxLI','PmaxLS')
%     fprintf('%5d %10.5f %10.5f %10.5f %10.5f %10.5f\n',ps.results.fv.varFolgaVW.')
%     fprintf('\n')
%     
%     fprintf('\n        --- Variáveis de fogla VV\n')
%     fprintf('%5s %10s \n', 'BarFV','xiC')
%     fprintf('%5d %10.5f \n',ps.results.fv.varFolgaVV.')
% end

% fprintf('\n----------- Fluxos de Potência\n');
% fprintf('   Perdas totais: %.4f MW e %.4f MVAr\n', ps.results.gerais(5), ps.results.gerais(6));
% fprintf('%5s %5s %5s %10s %10s %10s %10s %10s %10s\n','Nº','De','Para','Pde (MW)','Qde (MVAr)','Ppara','Qpara', 'Perda MW','Perda MVAr')
% fprintf('%5d %5d %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.ramo.')



end