%% Inicialización
Objeto = Tarea2('base_22.xls');
%% Causalidad A la Granger
[p,est_q] = Objeto.P1_RepAR(Objeto.raz_pi,15,10);

%%
Objeto.Causalidad_Granger(Objeto.raz_pi,1,Objeto.raz_y,5);