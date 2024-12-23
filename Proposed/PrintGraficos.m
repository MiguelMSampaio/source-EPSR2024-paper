%% ==== FUNC TO MAKE GRAPHS ===

function PrintGraficos(ps)

% % --- Voltage profile
% figure;
% hold on
% fA = plot(1:ps.Nb,ps.results.barra(:,2),'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
% set(gca, 'YGrid', 'on', 'XGrid', 'off')
% set(fA, 'MarkerFaceColor',get(fA,'color'))
% xlabel('Barra');
% ylabel('Tensão (pu)')
% set(gca, 'XTick',1:ps.Nb, 'XTickLabel',1:ps.Nb)
% 
% % -- Convergence curve
figure
semilogy(0:ps.results.gerais(1), ps.results.conv,'-x','LineWidth',1.5)
set(gca,'FontSize',10)
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
ylabel('$|f(x)|_{MAX}$','Interpreter','latex')
xlabel('$\upsilon$','Interpreter','latex')
% title('Erro x Iteração')
grid on

% if isfield(ps.results,'fv') && ps.results.fv.Nfv==1
%     figure
%     cmap = colormap(turbo(size(ps.results.resPrint,2)));
%     for k=1:size(ps.results.resPrint,2)
%         semilogy(0:ps.results.gerais(1), abs(ps.results.resPrint(:,k)),'Color', cmap(k, :))
%         hold on
%     end
%     hold off
%     legend('deltaG1','deltaG3')
%     if ps.fv.cont(1)==1
%         legend('deltaG1','deltaG3','deltaPhi1','deltaPhi2','deltaPhiFPF1','deltaPhiFPF2','deltaPhiFPF3')
%     elseif ps.fv.cont(1)==2
%         legend('deltaG1','deltaG3','deltaPhi2','deltaPhiVW1','deltaPhiVW2','deltaPhiVW3','deltaPhiVW4','deltaPhiVW5','deltaPhiVW6')
%     elseif ps.fv.cont(1)==3
%         legend('deltaG1','deltaG3','deltaPhi1','deltaPhi2','deltaPhiVV')
%     end
% 
%     figure
%     cmap = colormap(turbo(size(ps.results.varPrint,2)));
%     hold on
%     for k=1:size(ps.results.varPrint,2)
%         plot(0:ps.results.gerais(1), ps.results.varPrint(:,k),'Color',cmap(k,:))
%     end
%     hold off
%     legend('Vfv','Ifv','m','alpha')
%     if ps.fv.cont(1)==1
%         legend('Vfv','Ifv','m','alpha','xiC','Qli','Qls')
%     elseif ps.fv.cont(1)==2
%         legend('Vfv','Ifv','m','alpha','xiC','Pmax','Pls','PmaxLI','PmaxLS')
%     elseif ps.fv.cont(1)==3
%         legend('Vfv','Ifv','m','alpha','xiC')
%     end
% end

end