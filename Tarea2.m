classdef Tarea2
    %objeto desarrolado para la resoluci�n de la tarea 2
    %   
    
    properties
        df
        raz_pi
        raz_y
    end
    
    methods
        function obj = Tarea2(xlsx)
            %Importaci�n de datos de la tarea
            % df tiene 4 columnas, fecha, IMACEC, IPC y Tasa de pol�tica
            obj.df = readmatrix(xlsx);
            % Construir las series pi e y 
            % Se necesita perder las 12 primeras observaciones de cada
            % serie, MATLAB tiene el primer indice como 1, por tanto hay
            % que cortar 13:end
            div_pi = obj.df(13:end,3);
            div_y =  obj.df(13:end,2);
            % Luego de las �ltimas observaciones campoto tendr� del dato 12
            % meses adelante. Se deben perder 12 observaciones del vector a dividir.
            rez_pi = obj.df(1:end-12,3);
            rez_y = obj.df(1:end-12,2);
            %construyo la razon de las series por su �ltima observaci�n 12
            %meses atr�s
            raz_pi = div_pi ./ rez_pi;
            obj.raz_pi = log(raz_pi);
            raz_y = div_y ./ rez_y;
            obj.raz_y = log(raz_y);
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

