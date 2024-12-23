%% ========== DADOS DE ENTRADA DO SISTEMA ELÉTRICO

ps.sistema='IEEE 14 barras com FV';
ps.Sbase = 100e6;

%-------- Dados das barras (tipo: 1 - referência; 2 - PV; 3 - PQ; 4 - FV)
%OBS.: Considera-se que a barra de ref está na posição 1
%                                                         Geração                   Carga
%               Nº    Tipo      V(pu)       Ang(rad)    P(MW)     Q(MVAr)      P(MW)      Q(MVAr)      Bsh(pu)   Zona_Tensão  
%              (1)    (2)       (3)         (4)         (5)         (6)         (7)         (8)         (9)         (10)
ps.dadosBarra =[1       1       1.06        0           0           0           0           0           0           1;
%                 2       2       1.045       0           40          0           21.7        12.7        0           1;
                2       3       1.0         0           0           0           21.7        12.7        0           1;
%                 3       2       1.01        0           0           0           94.2        19          0           1;
                3       3       1.0         0           0           0           94.2        19          0           1;
                4       3       1           0           0           0           47.8        -3.9        0           1;
                5       3       1           0           0           0           7.6         1.6         0           1;
%                 6       2       1.07        0           0           0           11.2        7.5         0           2;
                6       3       1.0         0           0           0           11.2        7.5         0           2;
                7       3       1           0           0           0           0           0           0           0;
%                 8       2       1.09        0           0           0           0           0           0           3;
                8       3       1.0         0           0           0           0           0           0           3;
                9       3       1           0           0           0           29.5        16.6        0.19        2;
                10      3       1           0           0           0           9           5.8         0           2;
                11      3       1           0           0           0           3.5         1.8         0           2;
                12      3       1           0           0           0           6.1         1.6         0           2;
                13      3       1           0.00        0.00      0.00          13.50       5.80        0.00        1 ;  
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

fv.dadosPVSystem = [ ...  %(1) (2) (3) (4)   (5)  (6) (7)  (8)   (9)
                            2   40  44  .95   1   0  .98    .8    40;
                            3   40  44  1    1   0  .98     .9    48;
                            6   40  44  1    1   0  .98     1     55;
                            8   40  44  1    1   0  .98     .92   50];
                        
fv.PTcurve = [ ... %(barra)(T1,T2,T3,T4)  (P1,P2,P3,P4)  - em pu de Pnom
                    2       0,25,75,100    1.1125,1,.775,.55 ;
                    3       0,25,75,100    1.1125,1,.775,.55 ;
                    6       0,25,75,100    1.1125,1,.775,.55 ;
                    8       0,25,75,100    1.1125,1,.775,.55 ];   % contruido baseado no coef. de Pmax

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
                             2   -1   1  -1  0.92,0.98,1.02,1.08  1,0,0,-1     3    -1   1.08,1.1    1,0;
                             3   1   1  -1  0.92,0.98,1.02,1.08  1,0,0,-1     3    -1   1.08,1.1    1,0;
                             6   2   1  -1  0.92,0.98,1.02,1.08  1,0,0,-1     3    -1   1.08,1.1    1,0;
                             8   1   1  -1  0.92,0.98,1.02,1.08  1,0,0,-1     3    -1   1.08,1.1    1,0];
                         



            