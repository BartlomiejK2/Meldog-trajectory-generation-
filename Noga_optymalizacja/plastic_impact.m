%% Obliczanie impaktu z optymalizacji
function [f] = plastic_impact(q,dq,dq_,Fdt,l)
D = RANE_noga_2D(q,dq,zeros(3,1),zeros(3,1),zeros(3,1));
M1 = RANE_noga_2D(q,dq,[1;0;0],zeros(3,1),zeros(3,1))-D;
M2 = RANE_noga_2D(q,dq,[0;1;0],zeros(3,1),zeros(3,1))-D;
M3 = RANE_noga_2D(q,dq,[0;0;1],zeros(3,1),zeros(3,1))-D;
M = [M1,M2,M3];
J = Jacobi(q,l);
f = M*(dq-dq_) -J'*Fdt;
end

