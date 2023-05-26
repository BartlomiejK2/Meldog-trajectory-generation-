function [F_i_1_i, T_i_1_i, tau] = inverse_dynamics(m_i,J_com_i,r_0_i_11,r_0_i,r_com_i,ddr_com_i,R_0_i,w_0_i,dw_0_i,F_i_i_11,T_i_i_11,type,u)
%% Funkcja generująca siłe napędową (i), siły i momenty reakcji w parze kinamtycznej (i-1,i) na podstawie ruchu członu (i) i sił i momentów w parze (i,i+1)

% Przyspieszeni grawitacji
g = [0, -9.81, 0]';

% Siła w parze (i-1,i)
F_i_1_i = m_i*ddr_com_i - m_i*g + F_i_i_11;

% Moment w parze (i-1,i)
T_i_1_i = R_0_i*J_com_i*(transpose(R_0_i))*dw_0_i + skew_matrix(w_0_i)*R_0_i*J_com_i*(transpose(R_0_i))*w_0_i + T_i_i_11 +skew_matrix(r_0_i_11-r_com_i)*F_i_i_11 - skew_matrix(r_0_i-r_com_i)*F_i_1_i;

% Siła napędowa
if type == "p"
    tau = u'*F_i_1_i;
end
if type == "o"
    tau = u'*T_i_1_i;
end
end

