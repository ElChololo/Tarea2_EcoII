%% Inicialización
Objeto = Causalidad_Granger_2('base_22.xls');
%% Causalidad A la Granger
%% Generar representaciones AR que aseguren ruido blanco
[p_i,~, wht_ns_i] = Objeto.P1_RepAR(Objeto.tpm,10,20);
[p_y,~, wht_ns_y] = Objeto.P1_RepAR(Objeto.raz_y,10,20);
[p_pi,~, wht_ns_pi] = Objeto.P1_RepAR(Objeto.raz_pi,10,20);
%% Regresionar los residuos en rezagos para evaluar causalidades
% pi
[phi_sign_pi_y,orden_pi_y] = Objeto.Causalidad_Granger(wht_ns_pi,Objeto.raz_y,15)
[phi_sign_pi_i,orden_pi_i] = Objeto.Causalidad_Granger(wht_ns_pi,Objeto.tpm,15)
% y
[phi_sign_y_pi,orden_y_pi] = Objeto.Causalidad_Granger(wht_ns_y,Objeto.raz_pi,15)
[phi_sign_y_i,orden_y_i] = Objeto.Causalidad_Granger(wht_ns_y,Objeto.tpm,15)
% i
[phi_sign_i_pi,orden_i_pi] = Objeto.Causalidad_Granger(wht_ns_i,Objeto.raz_pi,15)
[phi_sign_i_y,orden_i_y] = Objeto.Causalidad_Granger(wht_ns_i,Objeto.raz_y,15)

