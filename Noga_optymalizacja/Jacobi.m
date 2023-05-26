%% Jakobian nogi
function J = Jacobi(q,l)
Om = [0, -1;
    1, 0];
l_1 = l(1);
l_2 = l(2);
alfa = q(2);
beta = q(3);
J = [[0;1],Om*rot(alfa)*[l_1;0]+Om*rot(alfa+beta)*[l_2;0],Om*rot(alfa+beta)*[l_2;0]];
end
function R = rot(a)
R = [cos(a), -sin(a);
    sin(a), cos(a)];
end
