%% Równania stanu dla symulacji z parą obrotową jako kontakt z podłożem
function [f] = dyn_simulation_noga_2D(q,dq,t,l)
C_G = RANE_noga_2D(q,dq,zeros(3,1),zeros(3,1),zeros(3,1));
M1 = RANE_noga_2D(q,dq,[1;0;0],zeros(3,1),zeros(3,1))-C_G;
M2 = RANE_noga_2D(q,dq,[0;1;0],zeros(3,1),zeros(3,1))-C_G;
M3 = RANE_noga_2D(q,dq,[0;0;1],zeros(3,1),zeros(3,1))-C_G;
M = [M1,M2,M3];
A = [M,-Jacobi(q,l)';
    Jacobi(q,l),zeros(2,2)];
b = [t-C_G;
    -dJacobi(q,dq,l)*dq];
x= A\b;
ddq = x(1:3,1);
f = [dq;ddq];
end

