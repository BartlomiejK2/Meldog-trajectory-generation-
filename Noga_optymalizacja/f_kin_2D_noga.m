function [xc,yc,dxc,dyc,ddxc,ddyc,Rc] = f_kin_2D_noga(q,dq,ddq)
l_1 = 0.25;
l_2 = 0.25;
xc = l_1*cos(q(2))+l_2*cos(q(2)+q(3));
yc = q(1)+l_1*sin(q(2))+l_2*sin(q(2)+q(3));
dxc = -l_1*sin(q(2))*dq(2)-l_2*sin(q(2)+q(3))*(dq(2)+dq(3));
dyc = dq(1)+l_1*cos(q(2))*dq(2)+l_2*cos(q(2)+q(3))*(dq(2)+dq(3));
ddxc = -l_1*sin(q(2))*ddq(2)-l_2*sin(q(2)+q(3))*(ddq(2)+ddq(3)) -l_1*cos(q(2))*(dq(2))^2-l_2*cos(q(2)+q(3))*(dq(2)+dq(3))^2;
ddyc = ddq(1) +l_1*cos(q(2))*ddq(2)+l_2*cos(q(2)+q(3))*(ddq(2)+ddq(3)) -l_1*sin(q(2))*(dq(2))^2-l_2*sin(q(2)+q(3))*(dq(2)+dq(3))^2;
Rc = rot(q(2))*rot(q(3));
end
function R = rot(a)
R = [cos(a), -sin(a);
    sin(a), cos(a)];
end