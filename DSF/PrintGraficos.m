%% ==== FUNC TO MAKE GRAPHS ===
% -- ENTRADA: ps [nested struct]  - power flow results
% -- SAIDA: print the graphs
function PrintGraficos(ps)

% --- Voltage profile
% figure;
% hold on
% fA = plot(1:ps.Nb,ps.results.barra(:,2),'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
% set(gca, 'YGrid', 'on', 'XGrid', 'off')
% set(fA, 'MarkerFaceColor',get(fA,'color'))
% xlabel('Barra');
% ylabel('Tensão (pu)')
% set(gca, 'XTick',1:ps.Nb, 'XTickLabel',1:ps.Nb)

% -- Convergence curve
% - inner loop error
figure
semilogy(1:ps.results.gerais(1), ps.results.conv,'-x','LineWidth',1.5)
set(gca,'FontSize',10)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
ylabel('$|g(x)|_{MAX}$','Interpreter','latex')
xlabel('Iter')
% title('Erro x Iteração')
grid on

% - Monitored voltage
% if ps.results.fv.NinvCtrl > 0
%     figure
%     plot(1:ps.results.gerais(7)+1,ps.results.fv.Vmon,'-x','LineWidth',1.5)
%     set(gca,'FontSize',10)
% %     title('Vmon x Iter Ctrl')
%     grid on
% end
    
end
