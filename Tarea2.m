classdef Tarea2
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
        
        function [p] = P1_RepAR(obj,serie,orden_max_ar,n_correlaciones)
            %% Estimar la mejor representación AR de la serie          
            for j=1:orden_max_ar
                % Perder observaciones(MODIFICABLE SI SE NECESITAN LOS AR CON MISMO NUMERO DE OBS)
                Var_Depe = serie(j:end);
                obs =length(Var_Depe);
                Regresores = obj.Get_Regresores(serie,j);
                % Regresión AR(orden)
                phis = (Regresores'*Regresores)\(Regresores'*serie);
                resid = serie - Regresores*phis;
                resid_2 = resid'*resid;
                % Computar test errores ruido blanco
                
                Test_pass = obj.Test_Q(resid,resid_2,n_correlaciones,obs);
                if Test_pass == 0
                    p=j;
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
        
        function Boolean_Test = Test_Q(obj,resid,resid_2,n_correlaciones,obs)
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
    end
end

