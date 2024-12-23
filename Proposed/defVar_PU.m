%% =============== STRUCTS DEFINITION AND TRANSFORMATION IN PU ======

% ------- bus data
Nb = length(ps.dadosBarra(:,1));                    % Quantidade de barras
DBAR.bar = ps.dadosBarra(:,1);                      % indice das barras
DBAR.tipo = ps.dadosBarra(:,2);                     % Tipo de cada barra (1 - ref, 2 - PV, 3 - PQ, 4 - fotovoltaico)
DBAR.V = ps.dadosBarra(:,3);                        % Tensão das barras em pu
DBAR.theta = ps.dadosBarra(:,4);                    % Ângulo das barras em radianos
DBAR.Pgen = 1e6*(ps.dadosBarra(:,5))/ps.Sbase;      % Potências de geração transformada em pu 
DBAR.Qgen = 1e6*(ps.dadosBarra(:,6))/ps.Sbase;
DBAR.Pfv = zeros(Nb,1);                             % Potências de geração de FV
DBAR.Qfv = zeros(Nb,1);
DBAR.Pc = 1e6*(ps.dadosBarra(:,7))/ps.Sbase;        % Potências de carga em pu
DBAR.Qc = 1e6*(ps.dadosBarra(:,8))/ps.Sbase;
DBAR.Pesp = DBAR.Pgen - DBAR.Pc;                    % Potências esperadas em cada barra em pu
DBAR.Qesp = DBAR.Qgen - DBAR.Qc;
DBAR.mhoSh = ps.dadosBarra(:,9);                    % Susceptância conectada unicamente a barra, em pu (ex.: banco cap)
DBAR.index_ref = find(DBAR.tipo==1);                % Endereço das barras de cada tipo
DBAR.index_PV = find(DBAR.tipo==2);        
DBAR.index_PQ = find(DBAR.tipo==3);    
DBAR.index_FV = find(DBAR.tipo==4);                         % fotovoltaico
DBAR.Nref = length(DBAR.index_ref);                 % Quantidade de barras de cada tipo
DBAR.NPV = length(DBAR.index_PV);
DBAR.NPQ = length(DBAR.index_PQ); 
DBAR.NFV = length(DBAR.index_FV);
DBAR.QespPQ = zeros(DBAR.NPQ,1);              %definição do vetor de pot reativo especificada apenas para as barras PQ
for k=1:DBAR.NPQ
    m=DBAR.index_PQ(k);
    DBAR.QespPQ(k) = DBAR.Qesp(m);
end


% -------- Line data
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


% -------- PV system data
if exist('fv','var')
    Nfv = length(fv.dados(:,1));            % Qntd de plantas FV
    % geral
    fv.bar = fv.dados(:,1);                 % indice da barra que o FV está
    fv.cont = fv.dados(:,2);                % 1: PQ/ 2: PV
    fv.index_FPF = find(fv.cont==1);        % Endereço das plantas FV com controle FPF
    fv.index_VW = find(fv.cont==2);         % " controle VW
    fv.index_VV = find(fv.cont==3);         %  " VV
    fv.Nfpf = length(fv.index_FPF);         % Número de FV com controle FPF
    fv.Nvw = length(fv.index_VW);           % " controle VW
    fv.Nvv = length(fv.index_VV);           %   " VV
    % inversor
    fv.Snom = fv.dados(:,3)/ps.Sbase;       % capacidade do inversor em pu [VA -> pu]
    fv.eta = fv.dados(:,4);                 % eficiencia do inversor
    % trafo acoplamento
    fv.r_t = fv.dados(:,5);                 % resistencia do trafo de acoplacamento
    fv.x_t = fv.dados(:,6);                 % reatancia " 
    fv.z_t = fv.r_t + 1i*fv.x_t;            % impedancia "
    fv.y_t = ones(Nfv,1)./fv.z_t;           % admitancia "
    fv.g_t = real(fv.y_t);                  % condutancia "
    fv.b_t = imag(fv.y_t);                  % susceptancia "
    % conf arranjo
    fv.Nss = fv.dados(:,7);                 % nº de modulos em serie
    fv.Npp = fv.dados(:,8);                % em paralelo
    % datasheet modulo no STC
    fv.Isc = fv.dados(:,9);                % corrente de curto circuito do modulo
    fv.Voc = fv.dados(:,10);                % tensao de circuito aberto
    fv.Vmpp = fv.dados(:,11);               % tensao MPP
    fv.Impp = fv.dados(:,12);               % corrente MPP
    fv.Pn = fv.Vmpp.*fv.Impp;                % capacidade do módulo FV
    fv.Ki = fv.dados(:,13);                 % coef. de temperatura de Isc
    fv.Kv = fv.dados(:,14);                 % " de Voc
    fv.ns = fv.dados(:,15);                 % nº de celular em serie no modulo
    % condições ambiente
    fv.Gin = fv.dados(:,16);                % irradiancia simulada, W/m2 (STC - 1000 W/m2)
    fv.Tin = fv.dados(:,17);                % temperatura simulada, ºC (STC - 25 ºC)
    
    % variaveis
    fv.Pkm = zeros(Nfv,1);                  % fluxo de potência através do trafo
    fv.Qkm = zeros(Nfv,1);
    fv.VbaseCC =  fv.Vmpp .* fv.Nss;        % tensão de base no elo CC 
    fv.Vbase = sqrt(3/8)*fv.VbaseCC;        % tensão de base na barra k (pós trafo)
%     fv.Vbase = ones(Nfv,1)*1e3;              % para o caso 5 barras e o teste IEEE 14
%     fv.Vbase = ones(Nfv,1)*516.65;      % Juarez IEEE 14 sem violação
%     fv.Vbase = ones(Nfv,1)*510.66;      % Juarez IEEE 14 violação
    
    % constante de idealidade do diodo
    for i=1:Nfv
        fv.a(i) = (fv.Vmpp(i)-fv.Voc(i))/( (fv.ns(i)*k_boltz*(fv.Tin(i)+273.15)/q) * log((fv.Isc(i)-fv.Impp(i))/fv.Isc(i)) );
    end

    
else
    Nfv = 0;
    fv = [];
end