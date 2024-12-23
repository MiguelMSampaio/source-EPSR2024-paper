%% ========== DADOS DE ENTRADA DO SISTEMA ELÉTRICO

ps.sistema='IEEE 14 barras com FV';
ps.Sbase = 100e6;

%-------- Dados das barras (tipo: 1 - referência; 2 - PV; 3 - PQ; 4 - FV)
%OBS.: Considera-se que a barra de ref está na posição 1
%                                                         Geração                   Carga
%               Nº    Tipo      V(pu)       Ang(rad)    P(MW)     Q(MVAr)      P(MW)      Q(MVAr)      Bsh(pu)   Zona_Tensão  
%              (1)    (2)       (3)         (4)         (5)         (6)         (7)         (8)         (9)         (10)
ps.dadosBarra =[1       1       1.06        0           0           0           0           0           0           1;
                2       2       1.045       0           40          0           21.7        12.7        0           1;
                3       2       1.01        0           0           0           94.2        19          0           1;
                4       3       1           0           0           0           47.8        -3.9        0           1;
                5       3       1           0           0           0           7.6         1.6         0           1;
                6       2       1.07        0           0           0           11.2        7.5         0           2;
                7       3       1           0           0           0           0           0           0           0;
                8       2       1.09        0           0           0           0           0           0           3;
                9       3       1           0           0           0           29.5        16.6        0.19        2;
                10      3       1           0           0           0           9           5.8         0           2;
                11      3       1           0           0           0           3.5         1.8         0           2;
                12      3       1           0           0           0           6.1         1.6         0           2;
                13      4       1           0           0           0           13.5        5.8         0           2;
                14      3       1           0           0           0           14.9        5           0           2];
                
% --------- Dados das linhas/ramos            
%                   Nº      De      Para    R           X           B/2         Tap-Setting (N1/N2)       :: tudo em pu;        
%                  (1)     (2)      (3)     (4)         (5)         (6)             (7)        %tap-setting é o a relação final de transformação do trafo N1/N2. definido pela posição do tap  
ps.dadosLinha = [   1       1       2       0.01938     0.05917     0.02640         1; 
                    2       1       5       0.05403     0.22304     0.02460         1;                
                    3       2       3       0.04699     0.19797     0.02190         1;             
                    4       2       4       0.05811     0.17632     0.01700         1;            
                    5       2       5       0.05695     0.17388     0.01730         1;
                    6       3       4       0.06701     0.17103     0.00640         1;
                    7       4       5       0.01335     0.04211     0.00000         1;
                    8       4       7       0.00000     0.20912     0.00000         0.978;
                    9       4       9       0.00000     0.55618     0.00000         0.969;
                    10      5       6       0.00000     0.25202     0.00000         0.932;
                    11      6       11      0.09498     0.19890     0.00000         1;
                    12      6       12      0.12291     0.25581     0.00000         1;
                    13      6       13      0.06615     0.13027     0.00000         1;
                    14      7       8       0.00000     0.17615     0.00000         1;
                    15      7       9       0.00000     0.11001     0.00000         1;
                    16      9       10      0.03181     0.08450     0.00000         1;
                    17      9       14      0.12711     0.27038     0.00000         1;
                    18      10      11      0.08205     0.19207     0.00000         1;
                    19      12      13      0.22092     0.19988     0.00000         1;
                    20      13      14      0.17093     0.34802     0.00000         1];    
            
%% ===== DADOS DO SISTEMA FV
% modo cont - 1: PQ / 2: PV / 3: VV
% Esp depende do modo cont 
%   - se PQ é FP esp; 
%   - se PV é Vref
% droop só pro PV. Se PQ, fica -1
% Nss, Npp - nº de módulos em série e paralelo do array
% Dados do datasheet do módulo (Isc,Voc,Vmpp,Impp,Ki,Kv,ns) em SI no STC
%   Ki, Kv em A/ºC e V/ºC, respectivamente
%   ns: nº de células em série do módulo
% Tin em ºC e Gin em W/m^2

% Dados do datasheet Kyocera KU330-8BCA 

%  Nºbarra  cont Snom[VA] eta_inv  R_T(pu) X_T(pu)  Nss   Npp   Isc   Voc    Vmpp  Impp   Ki    Kv    ns    Gin    Tin
%   (1)     (2)   (3)    (4)      (5)      (6)      (7)   (8)    (9)   (10)  (11)  (12)   (13)  (14)  (15)  (16)  (17)   

fv.dados = [ ...
    13       1       22e6      0.98      0       0.06    49   1236   8.74  50.3   40.7  8.11 .0052 -.1811  80    1000   25]; 


% - FPF - 1
%           Nºbarra     fator de potencia
fv.FPF = [   2               .8];

% - VW  - 2
% (P em pu de Snom)
%           Nºbarra     V1    P1    V2  P2
fv.VW = [  2            1.08   1     1.1  0];

% - VV - 3
% (Q em %Qdisp e no sentido de injeção)
%           Nºbarra     V1  Q1   V2  Q2   V3  Q3   V4  Q4
fv.VV = [  2          .92 1   .98  0    1.02 0   1.08 -1];



            