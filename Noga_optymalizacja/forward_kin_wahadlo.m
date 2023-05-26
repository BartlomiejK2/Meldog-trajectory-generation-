%% Kinematyka prosta dla końcówki nogi względem członu 1
function [r] = forward_kin_wahadlo(q,l)
l_1 = l(1);
l_2 = l(2);
x = l_1*cos(q(1))+l_2*cos(q(1)+q(2));
y = l_1*sin(q(1))+l_2*sin(q(1)+q(2));
r = [x;y];
end

