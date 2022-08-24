classdef Causalidad_Granger_2
    %objeto desarrolado para la resolución de la tarea 2
    %   
    
    properties
        df
        raz_pi
        raz_y
        tpm
    end
    
    methods
        function obj = Causalidad_Granger_2(xlsx)
            %Importación de datos de la tarea
            % df tiene 4 columnas, fecha, IMACEC, IPC y Tasa de política
            obj.df = readmatrix(xlsx);
            % serie tpm
            obj.tpm = obj.df(13:end,4);
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
        
        function [p,est_q,wht_ns] = P1_RepAR(obj,serie,orden_max_ar,n_correlaciones)
            %% Estimar la primera representación AR de la serie que genere ruidos blancos        
            for j=1:orden_max_ar
                % Perder observaciones para estimar el AR(orden_ar)
                Var_Depe = serie(j+1:end);
                obs =length(Var_Depe);
                Regresores = obj.Get_Regresores(serie,j);
                % Regresión AR(orden)
                phis = (Regresores'*Regresores)\(Regresores'*Var_Depe);
                resid = Var_Depe - Regresores*phis;
                resid_2 = resid'*resid;
                
                % Computar test errores ruido blanco
                
                [Test_pass,est_Q] = obj.Test_Q(resid,resid_2,n_correlaciones,obs);
                if Test_pass == 1
                    p=j;
                    est_q = est_Q ;
                    wht_ns = resid;
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
                vec_corr_est(k) = (vec_corr_est(k)^2) /  (obs-k);
            end
            est_Q = obs*(obs+2)*sum(vec_corr_est);
            % Devuelve 0 si rechazo la hipótesis nula, 1 si no lo rechazo
            Boolean_Test = est_Q < chi2inv(0.95,n_correlaciones);
        end
        
        function [phi_sig,orden] = Causalidad_Granger(obj,resid_wht,serie_2,orden_rez)
            % Ya obtenidos ruidos blancos, ver si los rezagos de la serie 2
            % tienen coeficientes asociados distintos de 0 en una regresión
            % contra esos ruidos blancos
            % Testear suficientes rezagos;
            ind = 0;
            for i=1:orden_rez
                Var_Depe = resid_wht;
                Regresores = obj.Get_Regresores(serie_2,i);
                % Dado que los errores vienen de procesos AR() se perderan
                % observaciones en la matriz de regresores. Además si el
                % orden de rezagos a utilizar es lo suficientemente grande
                % también se deberan perder observaciones en los
                % regresores.
                n_var_depe = length(Var_Depe);
                n_reg =length(Regresores);
                if n_var_depe < n_reg
                   diff = n_reg-n_var_depe;
                   
                   Regresores = Regresores(diff+1:end,:);
                elseif  n_var_depe > n_reg
                    diff = n_var_depe-n_reg;
                    Var_Depe=Var_Depe(diff+1:end);
                end
                coef_est = (Regresores'*Regresores)\(Regresores'*Var_Depe);
                resid = Var_Depe - Regresores*coef_est;
                resid_2 = resid'*resid;
                T =length(Var_Depe);
                k = length(coef_est);
                sigma_2_est = resid_2 / ( T-k );
                var_betas = sigma_2_est * ( (Regresores'*Regresores)^(-1) );
                % Computar la significancia estadística individual
                
                err_std = zeros(k,1);
                estad_t = zeros(k,1);
                valor_p = zeros(k,1);
                for j=2:k
                   err_std(j) = sqrt(var_betas(j,j));
                   t_est=coef_est(j)/err_std(j);
                   estad_t(j) = t_est;
                   valor_p(j) = tcdf(abs(t_est),T-k,'upper')*2;
                   if valor_p(j) < 0.05
                       %coeficiente significativo

                       
                       ind=1;
                       
                   end
                end
                phi_sig= [coef_est(2:end) err_std(2:end) estad_t(2:end) valor_p(2:end) ];
                if ind==1
                    fprintf("Significancia encontrada\n")
                    orden = i;
                    break
                end

            end
              
            
        end
    end
end
