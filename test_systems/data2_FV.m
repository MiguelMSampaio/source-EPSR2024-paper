
ps.sistema='2 barras exercício MONTICELLI com FV';

ps.Sbase = 100e6;           % 100 MVA - potencia de base

%-------- Dados das barras (tipo: 1 - referência; 2 - PV; 3 - PQ; 4 - FV)
%OBS.: Considera-se que a barra de ref está na posição 1
%                                                         Geração                   Carga
%               Nº    Tipo      V(pu)       Ang(rad)    P(MW)     Q(MVAr)      P(MW)      Q(MVAr)      Bsh(pu)   Zona_Tensão  
%              (1)    (2)       (3)         (4)         (5)         (6)         (7)         (8)         (9)         (10)
ps.dadosBarra =[1       1       0.95           0           0           0           0           0           0           1;
                2       4       1           0           0           0           0          0           0           1];

% --------- Dados das linhas/ramos            
%               Nº      De      Para    R           X           B/2         Tap-Setting         :: tudo em pu;        
%              (1)     (2)      (3)     (4)         (5)         (6)             (7)            %tap-setting é o a relação final de transformação do trafo N1/N2. definido pela posição do tap  
ps.dadosLinha = [  1     1       2       0         0.01           0             1];

%% ===== DADOS DO SISTEMA FV
% modo cont (Esp depende dele)
%  - 1: FPF
%  - 2: VW
%  - 3: VV
% Nss, Npp - nº de módulos em série e paralelo do array
% Dados do datasheet do módulo (Isc,Voc,Vmpp,Impp,Ki,Kv,ns) em SI no STC
%   Ki, Kv em A/K e V/K, respectivamente. Dá igual em ºC, por ser diferença
%   ns: nº de células em série do módulo
% Tin em ºC e Gin em W/m^2

% Dados do datasheet Kyocera KU330-8BCA e conf de 24 MWp e 25 MVA
%     2        3    25e6   0.98      0       0.9     24  1515*2  8.74  50.3  40.7   8.11 .0052 -.1811  80   1000  25];

%  Nºbarra  cont Snom[VA] eta_inv R_T(pu) X_T(pu)  Nss   Npp   Isc   Voc    Vmpp  Impp   Ki    Kv    ns    Gin    Tin
%   (1)     (2)   (3)    (4)      (5)      (6)      (7)  (8)    (9)   (10)  (11)  (12)   (13)  (14)  (15)  (16)  (17)  
fv.dados = [ ...        % Snom = 44 MVA e Pnom = 40 MWp
    2       3    2*22e6      0.98      0       0.06    49   2*1236   8.74  50.3   40.7  8.11 .0052 -.1811  80    1000   25];

% - FPF - 1
%           Nºbarra     fator de potencia
fv.FPF = [   2               1];
% - VW  - 2
% (P em pu de Snom)
%           Nºbarra     V1    P1    V2  P2
fv.VW = [  2            1.08   1     1.1  0];
% - VV - 3
% (Q em %Qdisp e no sentido de injeção)
%           Nºbarra     V1  Q1   V2  Q2   V3  Q3   V4  Q4
fv.VV = [  2          .92 1   .98  0    1.02 0   1.08 -1];



