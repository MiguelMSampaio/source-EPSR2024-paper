
ps.sistema='2 barras exercício MONTICELLI com FV';

ps.Sbase = 100e6;           % 100 MVA - potencia de base

%-------- Dados das barras (tipo: 1 - referência; 2 - PV; 3 - PQ; 4 - FV)
%OBS.: Considera-se que a barra de ref está na posição 1
%                                                         Geração                   Carga
%               Nº    Tipo      V(pu)       Ang(rad)    P(MW)     Q(MVAr)      P(MW)      Q(MVAr)      Bsh(pu)   Zona_Tensão  
%              (1)    (2)       (3)         (4)         (5)         (6)         (7)         (8)         (9)         (10)
ps.dadosBarra =[1       1       0.95        0           0           0           0           0           0           1;
                2       3       1           0           0           0           0          0           0           1];

% --------- Dados das linhas/ramos            
%               Nº      De      Para    R           X           B/2         Tap-Setting         :: tudo em pu;        
%              (1)     (2)      (3)     (4)         (5)         (6)             (7)            %tap-setting é o a relação final de transformação do trafo N1/N2. definido pela posição do tap  
ps.dadosLinha = [  1     1       2       0         0.01           0             1];

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
                            2   40  44  1    1   0  .98    1    25];
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



