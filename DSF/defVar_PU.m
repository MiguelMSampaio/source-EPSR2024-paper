%% =============== STRUCTS DEFINITION AND TRANSFORMATION IN PU ======

%-------- dos dados de barra
Nb = length(ps.dadosBarra(:,1));                    % Quantidade de barras
DBAR.bar = ps.dadosBarra(:,1);                      % indice das barras
DBAR.tipo = ps.dadosBarra(:,2);                     % Tipo de cada barra (1 - ref, 2 - PV, 3 - PQ, 4 - fotovoltaico)
DBAR.V = ps.dadosBarra(:,3);                        % Tensão das barras em pu
DBAR.theta = ps.dadosBarra(:,4);                    % Ângulo das barras em radianos
DBAR.Pgen = 1e6*(ps.dadosBarra(:,5))/ps.Sbase;      % Potências de geração transformada em pu 
DBAR.Qgen = 1e6*(ps.dadosBarra(:,6))/ps.Sbase;
DBAR.Pc = 1e6*(ps.dadosBarra(:,7))/ps.Sbase;        % Potências de carga em pu
DBAR.Qc = 1e6*(ps.dadosBarra(:,8))/ps.Sbase;
DBAR.Pesp = DBAR.Pgen - DBAR.Pc;                    % Potências esperadas em cada barra em pu
DBAR.Qesp = DBAR.Qgen - DBAR.Qc;
DBAR.mhoSh = ps.dadosBarra(:,9);                    % Susceptância conectada unicamente a barra, em pu (ex.: banco cap)
DBAR.index_ref = find(DBAR.tipo==1);                % Endereço das barras de cada tipo
DBAR.index_PV = find(DBAR.tipo==2);        
DBAR.index_PQ = find(DBAR.tipo==3);    
DBAR.Nref = length(DBAR.index_ref);                 % Quantidade de barras de cada tipo
DBAR.NPV = length(DBAR.index_PV);                     
DBAR.NPQ = length(DBAR.index_PQ); 
DBAR.QespPQ = zeros(DBAR.NPQ,1);              %definição do vetor de pot reativo especificada apenas para as barras PQ
for k=1:DBAR.NPQ
    m=DBAR.index_PQ(k);
    DBAR.QespPQ(k) = DBAR.Qesp(m);
end


% -------- dos dados de linha
Nramo = length(ps.dadosLinha(:,1));            % Quantidade de ramos
DRAM.i = ps.dadosLinha(:,1);                % Indices dos ramos
DRAM.de = ps.dadosLinha(:,2);                       % Indice da barra 'de'/'k' do ramo
DRAM.para = ps.dadosLinha(:,3);                     % Indice da barra 'para'/'m' do ramo
DRAM.z = ps.dadosLinha(:,4) + 1i*ps.dadosLinha(:,5);   % Impedância série em pu
DRAM.y = ones(Nramo,1)./DRAM.z;                          % Admitância série em pu
DRAM.g = real(DRAM.y);                                   % Condutância em pu
DRAM.b = imag(DRAM.y);                                   % Susceptância em pu
DRAM.Bsh = ps.dadosLinha(:,6);                      % Suceptância em derivação, já dividida por 2. Sem os elementos dos trafos
DRAM.a = ps.dadosLinha(:,7);                        % Tap do trafo (E1/E2  N1/N2)... 

if exist('fv','var')
    Nfv = size(fv.dadosPVSystem,1);                     % Numero de FV
    fv.bar = fv.dadosPVSystem(:,1);                     % Identificacao das barras do FV
    % - PVSystem
    fv.Pnom = fv.dadosPVSystem(:,2)/(ps.Sbase/1e6);     % Potencia nominal do arranjo FV, em pu
    fv.Snom = fv.dadosPVSystem(:,3)/(ps.Sbase/1e6);     % Potencia nominal do inversor, em pu
    fv.fp = fv.dadosPVSystem(:,4);              % Fator de potencia
    fv.wattPri = fv.dadosPVSystem(:,5);         % Prioridade de potencia (1 - reativa, 0 - ativa)
    fv.fpPri = fv.dadosPVSystem(:,6);           % Prioridade fator de potencia (tem prioridade sobre wattPri)
    fv.eta = fv.dadosPVSystem(:,7);             % Eficiencia do inversor
    fv.G = fv.dadosPVSystem(:,8);               % Irradiacao em kW/m2
    fv.T = fv.dadosPVSystem(:,9);               % Temperatura em C
    % - InvControl
    fv.cont = fv.dadosInvControl(:,2);              % Controle do FV (1 = Volt-Var - padrão / 2 = Volt-Watt / 3 = VV+VW / -1 = nao tem invcontrol)
    fv.temInvControl = ismember(fv.cont,[1 2 3]);       % veto Nfv,1 com 1 
    fv.NinvCtrl = nnz(fv.temInvControl);                % Numero de InvControls com VV, VW ou VV+VW
    fv.index_invCtrl = find(fv.temInvControl);          % enderecos das barras com InvCtrl
    fv.invControl = fv.dadosInvControl(:,3:end);    % Valores do invControl a partir de MODO. Se MODO=-1 eh controle FPF
    fv.FdeltaQ = fv.invControl([fv.index_invCtrl],2);   % Valor de entrada do deltaQ
    fv.FdeltaP = fv.invControl([fv.index_invCtrl],12);  % Valor de eentrada do deltaP
    for i=1:fv.NinvCtrl
        fvI = fv.index_invCtrl(i);
        if fv.invControl(fvI,2) == -1
            fv.invControl(fvI,2) = 0.5;
        end
        if fv.invControl(fvI,12) == -1
            fv.invControl(fvI,12) = 0.5;
        end
    end
    
    % - variaveis que precisam ser inicializadas
    fv.Pcc = zeros(Nfv,1);      % potencia ativa CC
    fv.Pca = zeros(Nfv,1);      % potencia ativa CA
    fv.Pca_lin = zeros(Nfv,1);  % potencia ativa CA desejada
    fv.Plimit = fv.Snom;        % padrao para a potencia ativa limite (do VW)
    fv;Plimit_min = zeros(Nfv,1);   % potencia ativa limite utiliza
    fv.Qca = zeros(Nfv,1);          % potencia reativa CA
    fv.Qca_lin = zeros(Nfv,1);      % potencia reativa desejada no CA
    fv.Qdesired = zeros(Nfv,1);     % potencia reativa desejada pelo InvControl (pos fator de escala)
    
    % invcontrol
    fv.Vmon = zeros(fv.NinvCtrl,1);     % tensao monitorada
    for i=1:fv.NinvCtrl
        fv.Vmon(i,1) = DBAR.V(fv.bar(fv.index_invCtrl(i)));
    end
    fv.Pbase = zeros(fv.NinvCtrl,1);    % potencia ativa de base (curva VW)
    fv.Pl_fun = zeros(fv.NinvCtrl,1);   % potencia ativa requiistada pela VW
    fv.Pl_lim = zeros(fv.NinvCtrl,1);   % potencia ativa limite pela capacidade do inversor
%     fv.Pl_end = zeros(fv.NinvCtrl,1);   % potencia ativa final (pre fator de escala)
    fv.Pl_end = fv.Plimit(fv.index_invCtrl);
%     fv.Pl_end = fv.Pnom(fv.index_invCtrl);
    fv.Qbase = zeros(fv.NinvCtrl,1);    % potencia reativa de base (curva VV)
    fv.Qd_fun = zeros(fv.NinvCtrl,1);   % potencia reativa requisitada pela VV
    fv.Qd_lim = zeros(fv.NinvCtrl,1);   % potencia reativa limite pela capacidade do inversor
    fv.Qd_end = zeros(fv.NinvCtrl,1);   % potencia reativa final (pre fator de escala)
    
end
                                            
