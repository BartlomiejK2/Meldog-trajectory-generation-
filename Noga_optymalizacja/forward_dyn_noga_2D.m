%% Dynamika w przód nogi
function [f] = forward_dyn_noga_2D(q,dq,t,Fz,Tz)
D = RANE_noga_2D(q,dq,zeros(3,1),Fz,Tz);
M1 = RANE_noga_2D(q,dq,[1;0;0],Fz,Tz)-D;
M2 = RANE_noga_2D(q,dq,[0;1;0],Fz,Tz)-D;
M3 = RANE_noga_2D(q,dq,[0;0;1],Fz,Tz)-D;
M = [M1,M2,M3];
ddq = M\(t-D);
f = [dq;ddq];
end

