%% Kinematyka prosta dla końcówki nogi względem układu 0
function [r] = forward_kin_noga_2D(q,l)
l_1 = l(1);
l_2 = l(2);
x = l_1*cos(q(2))+l_2*cos(q(2)+q(3));
y = q(1)+l_1*sin(q(2))+l_2*sin(q(2)+q(3));
r = [x;y];
end

