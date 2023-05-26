%% Jakobian nogi względem członu 1
function J = Jacobi_wahadlo(q,l)
Om = [0, -1;
    1, 0];
l_1 = l(1);
l_2 = l(2);
alfa = q(1);
beta = q(2);
J = [Om*rot(alfa)*[l_1;0]+Om*rot(alfa+beta)*[l_2;0],Om*rot(alfa+beta)*[l_2;0]];
end
function R = rot(a)
R = [cos(a), -sin(a);
    sin(a), cos(a)];
end
