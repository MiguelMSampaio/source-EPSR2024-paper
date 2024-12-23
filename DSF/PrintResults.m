%% ==== FUNC TO PRINT RESULTS ===== 
% -- ENTRADA: ps [nested struct] -> power flow results
% -- SAIDA: print on command window the results
function PrintResults(ps)

fprintf('  =================  Fluxo de Potência %s - Método Newton Raphson  =================== \n',ps.sistema);
fprintf('  Número de iterações      : %d \n', ps.results.gerais(1));
fprintf('  Número de it Ctrl        : %d \n', ps.results.gerais(7));
fprintf('  Resíduo máximo           : %g \n', ps.results.gerais(2));
% fprintf('  Tempo de execução total  : %.4f s\n', ps.results.gerais(3));
fprintf('  Tempo de execução método : %.4f s\n', ps.results.gerais(4));

% fprintf('\n----------- Variáveis de estado\n'); 
% fprintf('%5s %10s %10s %10s %10s %10s %10s\n','Barra','Tensão(pu)','Theta(deg)','Pgen(MW)','Qgen(MVAr)','Pc(MW)','Qc(MVAr)')
% fprintf('%5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.barra.')


% fprintf('\n----------- Fluxos de Potência\n');
% fprintf('   Perdas totais: %.4f MW e %.4f MVAr\n', ps.results.gerais(5), ps.results.gerais(6));
% fprintf('%5s %5s %5s %10s %10s %10s %10s %10s %10s\n','Nº','De','Para','Pde (MW)','Qde (MVAr)','Ppara','Qpara', 'Perda MW','Perda MVAr')
% fprintf('%5d %5d %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n',ps.results.ramo.')



end