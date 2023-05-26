%% Kinematyka odwrotna prędkości nogi względem układu 0
function [dq_a] = inv_vel_noga_2D(q,dy,l)
Om = [0, -1;
    1, 0];
l_1 = l(1);
l_2 = l(2);
alfa = q(2);
beta = q(3);
J1 = [0; 1];
J2 = Om*rot(alfa)*[l_1;0]+ Om*rot(alfa)*rot(beta)*[l_2;0];
J3 = rot(alfa)*Om*rot(beta)*[l_2;0];
J = [J2,J3];
dq_a = J\(-J1*dy);
end
function R = rot(a)
R = [cos(a), -sin(a);
    sin(a), cos(a)];
end