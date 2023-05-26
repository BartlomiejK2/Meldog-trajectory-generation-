%% Równania stanu dla symulacji, gdzie siła reakcji podłoża zamodelowana jest jako nieliniowa sprężyna
function [f] = dyn_simulation_noga_2D_spring(q,dq,t,l)
R = reaction_force_spring(q,dq,l);
D = RANE_noga_2D(q,dq,zeros(3,1),[R;0],zeros(3,1));
M1 = RANE_noga_2D(q,dq,[1;0;0],[R;0],zeros(3,1))-D;
M2 = RANE_noga_2D(q,dq,[0;1;0],[R;0],zeros(3,1))-D;
M3 = RANE_noga_2D(q,dq,[0;0;1],[R;0],zeros(3,1))-D;
M = [M1,M2,M3];
ddq = M\(t-D);
f = [dq;ddq];
end

