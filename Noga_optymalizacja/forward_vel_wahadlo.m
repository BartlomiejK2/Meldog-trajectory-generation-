%% Kinematyka prosta dla prędkości końcówki nogi względem członu 1
function [dr] = forward_vel_wahadlo(q,dq,l)
l_1 = l(1);
l_2 = l(2);
dx = -l_1*sin(q(1))*dq(1)-l_2*sin(q(1)+q(2))*(dq(1)+dq(2));
dy = l_1*cos(q(1))*dq(1)+l_2*cos(q(1)+q(2))*(dq(1)+dq(2));
dr = [dx;dy];
end

