

ps.sistema = "data30_FV.m";

ps.Sbase = 100e6;

ps.dadosBarra = [ 
       1        1    1.10     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
%        2        2    1.043     0.00    40.00     0.00    21.70    12.70     0.00        1 ; 
       2        3    1.000     0.00     0.00     0.00    21.70    12.70     0.00        1 ; 
       3        3    1.000     0.00     0.00     0.00     2.40     1.20     0.00        1 ; 
       4        3    1.000     0.00     0.00     0.00     7.60     1.60     0.00        1 ; 
       5        2    1.10     0.00     0.00     0.00    94.20    19.00     0.00        1 ; 
       6        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
       7        3    1.000     0.00     0.00     0.00    22.80    10.90     0.00        1 ; 
       8        2    1.10     0.00     0.00     0.00    30.00    30.00     0.00        1 ; 
       9        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      10        3    1.000     0.00     0.00     0.00     5.80     2.00     0.19        1 ; 
      11        2    1.082     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      12        3    1.000     0.00     0.00     0.00    11.20     7.50     0.00        1 ; 
      13        2    1.071     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      14        3    1.000     0.00     0.00     0.00     6.20     1.60     0.00        1 ; 
      15        3    1.000     0.00     0.00     0.00     8.20     2.50     0.00        1 ; 
      16        3    1.000     0.00     0.00     0.00     3.50     1.80     0.00        1 ; 
      17        3    1.000     0.00     0.00     0.00     9.00     5.80     0.00        1 ; 
      18        3    1.000     0.00     0.00     0.00     3.20     0.90     0.00        1 ; 
      19        3    1.000     0.00     0.00     0.00     9.50     3.40     0.00        1 ; 
      20        3    1.000     0.00     0.00     0.00     2.20     0.70     0.00        1 ; 
      21        3    1.000     0.00     0.00     0.00    17.50    11.20     0.00        1 ; 
      22        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      23        3    1.000     0.00     0.00     0.00     3.20     1.60     0.00        1 ; 
      24        3    1.000     0.00     0.00     0.00     8.70     6.70     0.04        1 ; 
      25        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      26        3    1.000     0.00     0.00     0.00     3.50     2.30     0.00        1 ; 
      27        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      28        3    1.000     0.00     0.00     0.00     0.00     0.00     0.00        1 ; 
      29        3    1.000     0.00     0.00     0.00     2.40     0.90     0.00        1 ; 
      30        3    1.000     0.00     0.00     0.00    10.60     1.90     0.00        1 ; 
 ];



% ===  DADOS DO SISTEMA FV (PVSystem e InvControl)

% -------- PVSystem
% 1 - ID -Identificação numérica da barra de conexão;
% 2 - P_NOM - Capacidade nominal do arranjo FV, em MW. Total, 'trifásico'
% 3 - S_NOM - Capacidade nominal do inversor (trifásica);
% 4 - FP - Fator de potência de saída. Padrão 1, caso especificado ativa controle fp fixo. Negativo (indutivo/absorve), positivo (capacitivo/fornece);
% 5 - WATTPriority - Prioridade potência ativa (1) ou reativa (0 - padrão);
% 6 - PFPriority - Prioriza o FP definido. Ligado (1), desligado (0). Tem prioridade sobre WATTPriority, exceto para controle VV;
% 7 - EFICIENCIA (eta) - Eficiência do inversor. Número entre 0 e 1;
% 8 - IRRADIANCIA (G) - Irradiância incidente no conjunto FV, em kW/m2;
% 9 - TEMPERATURA (T) - Temperatura do conjunto FV, na unidade de PTcurve.

fv.dadosPVSystem = [ ...  %(1) (2) (3) (4)      (5)  (6) (7)  (8)   (9)
                            2   40  44  0.8    1   0  .98    1    25];
fv.PTcurve = [ ... %(barra)(T1,T2,T3,T4)  (P1,P2,P3,P4)  - em pu de Pnom
                    2       0,25,75,100    1.1125,1,.775,.55];   % contruido baseado no coef. de Pmax

% ------- InvControl
% 1 - ID - Identificação numérica da barra q contém o PVSystem a ser controlado. TEM QUE TER UM PARA CADA Q TIVER
% 2 - MODO - Função inteligente a ser executada. Em código (1 = Volt-Var - padrão / 2 = Volt-Watt / 3 = VV+VW / -1 = nao tem invcontrol);
% 3 - REF_REACTIVE_POWER (refVAR) - Potência reativa base para o eixo Y da curva VV. (1 = VARAVAL - padrão / 2 = VARMAX - prop S_NOM do PVSystem);
% 4 - DELTA_Q - Limitador do passo da potência reativa. Entre 0 e 1;
% 5,6,7,8 - Vvv1,Vvv2,Vvv3,Vvv4 - Pontos no eixo X, de tensão, da curva VV;
% 9,10,11,12 - Qvv1,Qvv2,Qvv3,Qvv4 - Pontos no eixo Y da curva VV, EM PU;
% 13 - VOLTWATT_YAXIS (refWATT) - Potência ativa base para o eixo Y da curva VW. (1 = P_NOM - padrão / 2 = disponível / 3 = S_NOM);
% 14 - DELTA_P - Limitador do passo de potência ativa. Entre 0 e 1;
% 15,16 - Vvw1,Vvw2 - Pontos no eixo X, de tensão, da curva VW;
% 17,18 - Pvw1,Pvw2 - Pontos no eixo Y da curva VW, EM PU;

fv.dadosInvControl = [ ... %(1) (2) (3) (4) (5,6,7,8)           (9,10,11,12) (13) (14) (15,16)     (17,18)
                             2   1   1  -1  0.92,0.98,1.02,1.08  1,0,0,-1     3    -1   1.08,1.1    1,0];

% % ==== DADOS FV
% % Dados do datasheet Kyocera KU330-8BCA 
% %  Nºbarra  cont Snom[VA] eta_inv  R_T(pu) X_T(pu)  Nss   Npp   Isc   Voc    Vmpp  Impp   Ki    Kv    ns    Gin    Tin
% %   (1)     (2)   (3)    (4)      (5)      (6)      (7)   (8)    (9)   (10)  (11)  (12)   (13)  (14)  (15)  (16)  (17)   
% fv.dados = [ ...        % 20 MWp e 22 MVA  \\ 40 MWp e 44 MVA (2*)
%     2       1    2*22e6      0.98      0       0.06    49   2*1236   8.74  50.3   40.7  8.11 .0052 -.1811  80    1000   25];
% % --- FPF : 1
% %           Nºbarra     fator de potencia
% fv.FPF = [   2               0.95];
% % --- VW  : 2
% % (P em pu de Snom)
% %           Nºbarra     V1    P1    V2  P2
% fv.VW = [  2            1.08   1     1.1  0];
% % --- VV : 3
% % (Q em %Qdisp e no sentido de injeção)
% %           Nºbarra     V1  Q1   V2  Q2   V3  Q3   V4  Q4
% fv.VV = [  2          .92 1   .98  0    1.02 0   1.08 -1];


ps.dadosLinha = [ 
       1        1        2  0.01920  0.05750  0.02640  1.00000; 
       2        1        3  0.04520  0.16520  0.02040  1.00000; 
       3        2        4  0.05700  0.17370  0.01840  1.00000; 
       4        3        4  0.01320  0.03790  0.00420  1.00000; 
       5        2        5  0.04720  0.19830  0.02090  1.00000; 
       6        2        6  0.05810  0.17630  0.01870  1.00000; 
       7        4        6  0.01190  0.04140  0.00450  1.00000; 
       8        5        7  0.04600  0.11600  0.01020  1.00000; 
       9        6        7  0.02670  0.08200  0.00850  1.00000; 
      10        6        8  0.01200  0.04200  0.00450  1.00000; 
      11        6        9  0.00000  0.20800  0.00000  0.97800; 
      12        6       10  0.00000  0.55600  0.00000  0.96900; 
      13        9       11  0.00000  0.20800  0.00000  1.00000; 
      14        9       10  0.00000  0.11000  0.00000  1.00000; 
      15        4       12  0.00000  0.25600  0.00000  0.93200; 
      16       12       13  0.00000  0.14000  0.00000  1.00000; 
      17       12       14  0.12310  0.25590  0.00000  1.00000; 
      18       12       15  0.06620  0.13040  0.00000  1.00000; 
      19       12       16  0.09450  0.19870  0.00000  1.00000; 
      20       14       15  0.22100  0.19970  0.00000  1.00000; 
      21       16       17  0.05240  0.19230  0.00000  1.00000; 
      22       15       18  0.10730  0.21850  0.00000  1.00000; 
      23       18       19  0.06390  0.12920  0.00000  1.00000; 
      24       19       20  0.03400  0.06800  0.00000  1.00000; 
      25       10       20  0.09360  0.20900  0.00000  1.00000; 
      26       10       17  0.03240  0.08450  0.00000  1.00000; 
      27       10       21  0.03480  0.07490  0.00000  1.00000; 
      28       10       22  0.07270  0.14990  0.00000  1.00000; 
      29       21       22  0.01160  0.02360  0.00000  1.00000; 
      30       15       23  0.10000  0.20200  0.00000  1.00000; 
      31       22       24  0.11500  0.17900  0.00000  1.00000; 
      32       23       24  0.13200  0.27000  0.00000  1.00000; 
      33       24       25  0.18850  0.32920  0.00000  1.00000; 
      34       25       26  0.25440  0.38000  0.00000  1.00000; 
      35       25       27  0.10930  0.20870  0.00000  1.00000; 
      36       28       27  0.00000  0.39600  0.00000  0.96800; 
      37       27       29  0.21980  0.41530  0.00000  1.00000; 
      38       27       30  0.32020  0.60270  0.00000  1.00000; 
      39       29       30  0.23990  0.45330  0.00000  1.00000; 
      40        8       28  0.06360  0.20000  0.02140  1.00000; 
      41        6       28  0.01690  0.05990  0.00650  1.00000; 
];