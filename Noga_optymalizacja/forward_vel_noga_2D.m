%% Kinematyka prosta dla prędkości końcówki nogi względem układu 0
function [dr] = forward_vel_noga_2D(q,dq,l)
l_1 = l(1);
l_2 = l(2);
dx = -l_1*sin(q(2))*dq(2)-l_2*sin(q(2)+q(3))*(dq(2)+dq(3));
dy = dq(1)+l_1*cos(q(2))*dq(2)+l_2*cos(q(2)+q(3))*(dq(2)+dq(3));
dr = [dx;dy];
end

