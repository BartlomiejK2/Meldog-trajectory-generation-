function [tau] = RANE_noga_2D(q,dq,ddq,Fz,Tz)

%% Funkcja do generowania momentów napędowych dla nogi 2D
y = q(1);
alfa = q(2);
beta = q(3);
dy = dq(1);
dalfa = dq(2);
dbeta = dq(3);
ddy = ddq(1);
ddalfa = ddq(2);
ddbeta = ddq(3);


%% Dane:

% Człon 1:
m_1 = 2;
h = 0.2;
d = 0.1;
w = 0.2;
J_1 = [1/12*m_1*(h^2+d^2) 0 0;
    0 1/12*m_1*(d^2+w^2) 0;
    0 0 1/12*m_1*(h^2+d^2)];
r_com1_const = [0;0;0];
% Człon 2:
m_2 = 0.5;
l_1 = 0.25;
r_1 = 0.01;
r_com2_const = [l_1/2;0;0];
j_2 = 1/12*m_2*((l_1)^2+3*r_1^2);
J_2 = [1/2*m_2*r_1^2 0 0;
    0 j_2 0;
    0 0 j_2];

% Człon 3:
m_3 = 0.5;
l_2 = 0.25;
r_2 = 0.01;
r_com3_const = [l_2/2;0;0];
j_3 = 1/12*m_2*((l_2)^2+3*r_2^2);
J_3 = [1/2*m_2*r_2^2 0 0;
    0 j_3 0;
    0 0 j_3];

% Wektory między układami współrzędnych:
r_0_1_const = [0;0;0];
r_1_2_const = [0;0;0];
r_2_3_const = [l_1;0;0];
r_3_4_const = [l_2;0;0];

%% Kinematyka w przód:

[r_0_1, dr_0_1 ,ddr_0_1 ,R_0_1, w_0_1, dw_0_1, r_com_1, dr_com_1, ddr_com_1, u1] = forward_kinematics_pos(zeros(3,1),zeros(3,1),zeros(3,1),eye(3,3),zeros(3,1),zeros(3,1),r_com1_const,"y",y,dy,ddy,r_0_1_const,eye(3,3));
[r_0_2, dr_0_2 ,ddr_0_2 ,R_0_2, w_0_2, dw_0_2, r_com_2, dr_com_2, ddr_com_2, u2] = forward_kinematics_obr(r_0_1,dr_0_1,ddr_0_1,R_0_1,w_0_1,dw_0_1,r_com2_const,"z",alfa,dalfa,ddalfa,r_1_2_const,eye(3,3));
[r_0_3, dr_0_3 ,ddr_0_3 ,R_0_3, w_0_3, dw_0_3, r_com_3, dr_com_3, ddr_com_3, u3] = forward_kinematics_obr(r_0_2,dr_0_2,ddr_0_2,R_0_2,w_0_2,dw_0_2,r_com3_const,"z",beta,dbeta,ddbeta,r_2_3_const,eye(3,3));
[r_0_4, dr_0_4 ,ddr_0_4 ,R_0_4, w_0_4, dw_0_4, r_com_4, dr_com_4, ddr_com_4, u4] = forward_kinematics_obr(r_0_3,dr_0_3,ddr_0_3,R_0_3,w_0_3,dw_0_3,zeros(3,1),"z",0,0,0,r_3_4_const,eye(3,3));

%% Dynamika odwrotna:
[F_2_3, T_2_3, tau3] = inverse_dynamics(m_3,J_3,r_0_4,r_0_3,r_com_3,ddr_com_3,R_0_3,w_0_3,dw_0_3,-Fz,-Tz,"o",u3);
[F_1_2, T_1_2, tau2] = inverse_dynamics(m_2,J_2,r_0_3,r_0_2,r_com_2,ddr_com_2,R_0_2,w_0_2,dw_0_2,F_2_3,T_2_3,"o",u2);
[F_0_1, T_0_1, tau1] = inverse_dynamics(m_1,J_1,r_0_2,r_0_1,r_com_1,ddr_com_1,R_0_1,w_0_1,dw_0_1,F_1_2,T_1_2,"p",u1);


%% Siły napędowe:

tau = [tau1;tau2;tau3];

end

