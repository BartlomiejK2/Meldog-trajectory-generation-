%% Funkcja obliczająca siłe reakcji na podstawie równania ruchu
function [f] = calc_F(q,dq,ddq,t,l)
C_G = RANE_noga_2D(q,dq,zeros(3,1),zeros(3,1),zeros(3,1));
M1 = RANE_noga_2D(q,dq,[1;0;0],zeros(3,1),zeros(3,1))-C_G;
M2 = RANE_noga_2D(q,dq,[0;1;0],zeros(3,1),zeros(3,1))-C_G;
M3 = RANE_noga_2D(q,dq,[0;0;1],zeros(3,1),zeros(3,1))-C_G;
M = [M1,M2,M3];
F = (Jacobi(q,l)')\(M*ddq+C_G-t);
f = F;
end


