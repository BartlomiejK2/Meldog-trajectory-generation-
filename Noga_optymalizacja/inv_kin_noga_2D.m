%% Kinematyka odwrotna nogi względem układu 0
function [q_a] = inv_kin_noga_2D(y,p,l)
l_1 = l(1);
l_2 = l(2);
beta = -acos((y^2+p^2-l_1^2-l_2^2)/(2*l_1*l_2));
alfa_y = (-p*(l_2*sin(beta))-y*(l_1+l_2*cos(beta)))/(l_1^2+l_2^2+2*l_1*l_2*cos(beta));
alfa_x = p/(l_1+l_2*cos(beta))+l_2*sin(beta)/(l_1+l_2*cos(beta))*alfa_y;
alfa = atan2(alfa_y,alfa_x);
if(beta>0)
    beta = acos((y^2+p^2-l_1^2-l_2^2)/(2*l_1*l_2));
    alfa_y = (-p*(l_2*sin(beta))-y*(l_1+l_2*cos(beta)))/(l_1^2+l_2^2+2*l_1*l_2*cos(beta));
    alfa_x = p/(l_1+l_2*cos(beta))+l_2*sin(beta)/(l_1+l_2*cos(beta))*alfa_y;
    alfa = atan2(alfa_y,alfa_x);
end
q_a(1,1) = alfa;
q_a(2,1) = beta;
end

