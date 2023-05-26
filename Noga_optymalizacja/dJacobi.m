%% Funkcja generująca pochodną względem czasu Jakobianu nogi 
function dJ = dJacobi(q,dq,l)
Om = [0, -1;
    1, 0];
l_1 = l(1);
l_2 = l(2);
alfa = q(2);
beta = q(3);
dalfa = dq(2);
dbeta = dq(3);
dJ = [[0;0],Om*Om*rot(alfa)*[l_1;0]+Om*Om*rot(alfa+beta)*[l_2;0],Om*Om*rot(alfa+beta)*[l_2;0]]*dalfa +[[0;0],Om*Om*rot(alfa+beta)*[l_2;0],Om*Om*rot(alfa+beta)*[l_2;0]]*dbeta;
end
function R = rot(a)
R = [cos(a), -sin(a);
    sin(a), cos(a)];
end
