function [r_0_i, dr_0_i ,ddr_0_i ,R_0_i, w_0_i, dw_0_i, r_com_i, dr_com_i, ddr_com_i, u] = forward_kinematics(r_0_i_1,dr_0_i_1,ddr_0_i_1,R_0_i_1,w_0_i_1,dw_0_i_1,r_i_1_i,dr_i_1_i,ddr_i_1_i,R_i_1_i,w_i_1_i,dw_i_1_i,s_com, w)
%% Funkcja generująca położenie i orientacje układu wspołrzędnych i oraz położenie, prędkość i przyspieszenie środka masy członu i na podstawie danych min. z układu i-1 (tutaj i_1)

% Położenie układu
r_0_i = r_0_i_1 + R_0_i_1 * r_i_1_i;

% Prędkość układu
dr_0_i = dr_0_i_1 + skew_matrix(w_0_i_1)*R_0_i_1 * r_i_1_i + R_0_i_1 * dr_i_1_i;

% Przyspieszenie układu
ddr_0_i = ddr_0_i_1 + skew_matrix(dw_0_i_1)*R_0_i_1 * r_i_1_i + skew_matrix(w_0_i_1)*skew_matrix(w_0_i_1)*R_0_i_1 * r_i_1_i + 2*skew_matrix(w_0_i_1)*R_0_i_1 * dr_i_1_i+ R_0_i_1 * ddr_i_1_i;

% Orientacja układu
R_0_i = R_0_i_1*R_i_1_i;

% Prędkość kątowa układu
w_0_i = w_0_i_1 + R_0_i_1*w_i_1_i;

% Przyspieszenie kątowe układu
dw_0_i = dw_0_i_1 +skew_matrix(dw_0_i_1)*R_0_i_1*w_i_1_i + R_0_i_1*dw_i_1_i;

% Położenie CoM
r_com_i = r_0_i + R_0_i * s_com;

% Prędkość CoM
dr_com_i = dr_0_i + skew_matrix(w_0_i)*R_0_i*s_com;

% Przyspieszenie CoM

ddr_com_i = ddr_0_i + skew_matrix(dw_0_i)*R_0_i*s_com + skew_matrix(w_0_i)*skew_matrix(w_0_i)*R_0_i*s_com;

% Wersor pary kinematycznej

u = R_0_i*w;

end


