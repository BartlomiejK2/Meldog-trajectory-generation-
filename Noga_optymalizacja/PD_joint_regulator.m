%% Funkcja generująca moment napędowy z regulatora złączowego PD
function [t] = PD_joint_regulator(q_r,q,dq_r,dq) % q_r, dq_r - współrzędne i prędkości złączowe referencyjne, q, dq - współrzędne i prędkości złączowe rzeczywiste
K_1 = 4; %4.5;
K_2 = 4; %4.5;
D_1 = 0.1; %1;
D_2 = 0.1; %1;
t1 = K_1*(q_r(1)-q(1))+D_1*(dq_r(1)-dq(1));
t2 = K_2*(q_r(2)-q(2))+D_2*(dq_r(2)-dq(2));
t = [t1;t2];
end

