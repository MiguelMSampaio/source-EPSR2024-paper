%% ====== CLASS SPLINE SMOOTHING
classdef SplineSuav < handle
    
%     Métodos
%           SplineSuav - definição. define o tipo de spline e seus parâmetros
%           output - saida para Vm, a depender do tipo
%           del - derivada para Vm, a depender do tipo

    properties
        % -- gerais
        Vscale;         % escala de V
        Qscale;         % escala de Q
        tipo;           % tipo de spline ('cubic' ou 'circular')
        v1; v2;         % tensoes que inicia e finaliza a suavizacao
        
        % -- circular
        xc;         
        yc;
        r2;
        sign;
        
        % -- cubic
        xs;
        A;
        B;
        C;
        D;
        
    end
   
    methods
        %% === CONSTUTOR
        function e = SplineSuav(va,v0,vb, qa,q0,qb, ps)         % constroi o objeto
            % PARAMETROS:       va,v0,vb - valores de V
            %                   qa,q0,qb - valores de Q
            %                   ps - dados do sistema, contem Sbase e tolerancia
            % SAIDA:            e - objeto SplineSuav criado. Com os parametros de
            %                       circular ou cubic
            
            if nargin
            dva = v0-va;
            dvb = vb-v0;
            dqa = q0-qa;
            dqb = qb-q0;
            
            if (abs(dva) < abs(dvb)) 
                e.Vscale = 0.10*abs(dva);
            else
                e.Vscale = 0.10*abs(dvb);
            end
            if e.Vscale < 0.001 
                e.Vscale = 0.001;
            end
            
            e.Qscale = 0.10*( abs(dqa) + abs(dqb) );
            if e.Qscale < ps.tol 
                e.Qscale = ps.tol;
            end
            
            dya = dqa/e.Qscale;
            dyb = dqb/e.Qscale;
            dxa = dva/e.Vscale;
            dxb = dvb/e.Vscale;
            da = sqrt( dxa^2 + dya^2 );
            db = sqrt( dxb^2 + dyb^2 );
            
            % We need the length in the transformed coordinates to be > 2.0
            % otherwise point 1 or 2 will be more than 50% down segment
            if da < db 
                dMin = da;
            else
                dMin = db;
            end
            if dMin < 2 
                e.Qscale = e.Qscale * dMin/2;
                e.Vscale = e.Vscale * dMin/2;
                da = da * 2/dMin;
                db = db * 2/dMin;
                dxa = dxa * 2/dMin;
                dya = dya * 2/dMin;
                dxb = dxb * 2/dMin;
                dyb = dyb * 2/dMin;
            end
            
            percA = 1/da;
            percB = 1/db;
            x0 = v0/e.Vscale;
            y0 = q0/e.Qscale;
            x1 = x0 - dxa*percA;
            y1 = y0 - dya*percA;
            x2 = x0 + dxb*percB;
            y2 = y0 + dyb*percB;
            
            e.v1 = x1*e.Vscale;
            e.v2 = x2*e.Vscale;
            
            angleA = atan2(dya, dxa);
            angleB = atan2(dyb, dxb);
            if abs(angleB - angleA) > pi/15 % at 12 degrees, switch
                e.tipo = 'circular';
            else
                e.tipo = 'cubic';
            end
            
            % ---- Definicao dos parametros do spline  
            if strcmp(e.tipo,'circular')
%                 fprintf('circular\n')
                if y0 == y1
                    e.xc = x1;
                    e.yc = y2 - ((x0-x2)/(y0-y2))*(x1-x2);
                    e.r2 = (y1-e.yc)^2;
                    e.sign = +1;
                elseif y0 == y2 
                    e.xc = x2;
                    e.yc = y1 - ((x0-x1)/(y0-y1))*(x2-x1);
                    e.r2 = (y2-e.yc)^2;
                    e.sign = -1;
                elseif x0 == x1
                    e.yc = y1;
                    e.xc = x2 - ((y0-y2)/(x0-x2))*(y1-y2);
                    e.r2 = (x1-e.xc)^2;
                    e.sign = -1;
                elseif x0 == x2
                    e.yc = y2;
                    e.xc = x1 - ((y0-y1)/(x0-x1))*(y2-y1);
                    e.r2 = (x2-e.xc)^2;
                    e.sign = +1;
                else 
                    m1 = (y0-y1)/(x0-x1);
                    m2 = (y0-y2)/(x0-x2);
                    e.xc = (m1*(m2*y2 + x2) - m2*(m1*y1 + x1))/(m1-m2);
                    e.yc = ((m1*y1 + x1) - (m2*y2 + x2))/(m1-m2);
                    e.r2 = (x1-xc)^2 + (y1-yc)^2;
                    if m1 < m2 
                        e.sign = -1;
                    else
                        e.sign = +1;
                    end
                end
            elseif strcmp(e.tipo,'cubic')
                e.xs = x1;
                x0 = x0 - e.xs;
                x1 = x1 - e.xs;
                x2 = x2 - e.xs;
                % calculation assumes that the x values are increasing x1 < x0 < x2
                m1 = (y0-y1)/(x0-x1);
                m2 = (y0-y2)/(x0-x2);
                % Calculate values on shifted coordinates around [(x1+x2)/2, y1]
                at = (y1+y2)/2 - ((x1-x2)*(m1-m2))/8;
                bt = (3*(y1-y2))/(2*(x1-x2)) - (m1+m2)/4;
                ct = (m1-m2)/(2*(x1-x2));
                dt = (m1+m2)/((x1-x2)^2) - (2*(y1-y2))/((x1-x2)^3);
                % shift back to the coordinates of centered around original x1, y1
                x12 = (x1+x2)/2;
                e.A = at - bt*x12 + ct*x12^2 - dt*x12^3;
                e.B = bt - 2*ct*x12 + 3*dt*x12^2;
                e.C = ct - 3*dt*x12;
                e.D = dt;
            end
        
            end
        end
        
        %% ==== CALCULA O VALOR DA SPLINE
        function result = output(spline, Vm)
            % PARAMETROS:       spline - objeto SplineSuav
            %                   Vm - nivel de tensao 
            % SAIDA:            result - resultado do spline em Vm
%             fprintf('%.3f\n',spline.v1)
            if strcmp(spline.tipo,'circular')
                x = Vm/spline.Vscale;
                result = spline.r2 - (x - spline.xc)^2;
                if result > 0 
                    result = spline.yc + spline.sign * sqrt( result );
                else        % function should never be used in this range!
                    result = spline.yc;
                end
                result = result*spline.Qscale;
                
            elseif strcmp(spline.tipo,'cubic')
                x = Vm/spline.Vscale - spline.xs;
                result = spline.A + spline.B*x + spline.C*x^2 + spline.D*x^3;
                result = spline.Qscale*result;
            end
            
        end
        
        %% ====== CALCULA O VALOR DA DERIVADA EM Vm
        function result = del(spline,Vm)
            % PARAMETROS:       spline - objeto SplineSuav com a suavizacao
            %                   Vm - tensao a ser calculada da curva VV
            % SAIDA:            result - valor da derivada em Vm
%             fprintf('%.3f\n',spline.v1)
            if strcmp(spline.tipo,'circular')
                x = Vm/spline.Vscale;
                result = spline.r2 - (x - spline.xc)^2;
                if result > 0
                    result = -spline.sign*(x - spline.xc)/sqrt( result );
                else
                    result = -1e6;
                end
                result = result*(spline.Qscale/spline.Vscale);
            elseif strcmp(spline.tipo,'cubic')
                x = Vm/spline.Vscale - spline.xs;
                result = spline.B + 2*spline.C*x + 3*spline.D*x^2;
                result = result*(spline.Qscale/spline.Vscale);
            end
        end
        
       
    end
end