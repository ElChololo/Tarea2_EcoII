classdef Causalidad_Granger
    %objeto desarrolado para la resolución de la tarea 2
    %   
    
    properties
        df
        raz_pi
        raz_y
    end
    
    methods
        function obj = Tarea2(xlsx)
            %Importación de datos de la tarea
            % df tiene 4 columnas, fecha, IMACEC, IPC y Tasa de política
            obj.df = readmatrix(xlsx);
            % Construir las series pi e y 
            % Se necesita perder las 12 primeras observaciones de cada
            % serie, MATLAB tiene el primer indice como 1, por tanto hay
            % que cortar 13:end
            div_pi = obj.df(13:end,3);
            div_y =  obj.df(13:end,2);
            % Luego de las últimas observaciones campoto tendré del dato 12
            % meses adelante. Se deben perder 12 observaciones del vector a dividir.
            rez_pi = obj.df(1:end-12,3);
            rez_y = obj.df(1:end-12,2);
            %construyo la razon de las series por su última observación 12
            %meses atrás
            raz_pi = div_pi ./ rez_pi;
            obj.raz_pi = log(raz_pi);
            raz_y = div_y ./ rez_y;
            obj.raz_y = log(raz_y);
        end
        
        function [p,est_q] = P1_RepAR(obj,serie,orden_max_ar,n_correlaciones)
            %% Estimar la mejor representación AR de la serie          
            for j=1:orden_max_ar
                % Perder observaciones(MODIFICABLE SI SE NECESITAN LOS AR CON MISMO NUMERO DE OBS)
                Var_Depe = serie(j+1:end);
                obs =length(Var_Depe);
                Regresores = obj.Get_Regresores(serie,j);
                % Regresión AR(orden)
                phis = (Regresores'*Regresores)\(Regresores'*Var_Depe);
                resid = Var_Depe - Regresores*phis;
                resid_2 = resid'*resid;
                % Computar test errores ruido blanco
                
                [Test_pass,est_Q] = obj.Test_Q(resid,resid_2,n_correlaciones,obs);
                if Test_pass == 0
                    p=j;
                    est_q = est_Q ;
                    break
                end
            end
        end
        
        function Regresores = Get_Regresores(obj,serie,orden)
            %% Objeto para crear regresores de estimaciones AR,incluyendo constante
            obs = length(serie)-orden;
            Regresores= zeros(obs,orden+1);
            cte = ones(obs,1);
            Regresores(:,1) = cte;
            for i=1:orden
               rezagos = circshift(serie,i);
               Regresores(:,i+1) = rezagos(orden+1:end);
            end
        end
        
        function [Boolean_Test,est_Q] = Test_Q(obj,resid,resid_2,n_correlaciones,obs)
            %% Computar las correlaciones estimadas
            % Armar matriz que albergara las m primeras correlaciones
            vec_corr_est = zeros(n_correlaciones,1);
            for j=1:n_correlaciones
               vec_corr_est(j) =  (resid(j+1:end)' * resid(1:end-j)) / resid_2;
            end
            %Armar el estadístico Q
            for k=1:n_correlaciones
                vec_corr_est(j) = (vec_corr_est(j)^2) /  (obs-n_correlaciones);
            end
            est_Q = obs*(obs+2)*sum(vec_corr_est);
            Boolean_Test = est_Q < chi2inv(0.975,2);
        end
        
        function coef_est = Causalidad_Granger(obj,serie,p,serie_2,orden_rez)
            %Primero Armar el AR que hace que los errores sean ruido blanco
            Var_Depe= serie(p+1:end);
            Regresores_Ar = obj.Get_Regresores(serie,p);
            % Agregar los regresores externos
            for j=1:orden_rez
                Regresores_Ext = obj.Get_Regresores(serie_2,j);
                % Concatenar regresores, ojo que los regresores ext tiene
                % intercepto
                Regresores_Ext = Regresores_Ext(:,2:end);
                Regresores = [Regresores_Ar Regresores_Ext];
                %coef_est es un vector columna, con las 2 primeras filas
                %correspondientes al intercepto y al rezago 1
                coef_est = (Regresores'*Regresores)\(Regresores'*Var_Depe);
                % Realizar test para ver si los regresores Ext -> Externos
                % no son estadisticamente diferentes de 0.
                resid = Var_Depe - Regresores*coef_est;
                resid_2 = resid'*resid;
                T =length(Var_Depe);
                k = length(coef_Est);
                sigma_est = resid_2 / ( T-k );
                var_betas = sigma_est * ( (Regresores'*Regresores)^(-1) );
                if j ==1
                    est_t = coef_est/ var_betas(3,3);
                    val_critico =  tinv(0.95,(T - k);
                    if est_t > val_critico
                       %se rechaza la hipotesis nula que serie 2 no causa a serie 1 en primer rezago 
                    end
                else
                    %Armar el test F
                    
                    R = ones(k,1);
                    %% Ojo que estas asumiendo aqui que la  mejor representacion es un AR(1) para la serie 1
                    R(1:2) = [0;0];
                    c = zeros(k-2);
                    aux_f = (R'*coef_Est-c);
                    %Nuevamente ese k-2 asume que la mejor representación
                    %es AR(1)
                    est_f = ( (aux_f)'*[sigma_est*R'*(Regresores'*Regresores)/R]*aux_f ) / (k-2);
                    valor_critico = finv(0.95,k-2,T-k);
                    if est_f > valor_critico
                        %se rechaza la hipotesis nula
                        
                    end
                end
                
                
                
            end
            
        end
    end
end

