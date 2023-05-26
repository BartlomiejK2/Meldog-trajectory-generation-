%% Funkcja generująca siłę reakcji na podstawie modelu nieliniowej sprężyny i tłumika
function [R] = reaction_force_spring(q,dq,l)
% Siła Fy:
c_max = 50;
h = 1e-5;
k  =5*1e+4;
p = 1.1; % dla gumy
% Siła Fx:
u_s = 0.7;
u_d = 0.4;
v_s = 1e-1;
v_d = 1e-0;
v = Jacobi(q,l)*dq;
dx_k = v(1);
dy_k = v(2);
r =  forward_kin_noga_2D(q,l);
y_k = r(2);
if(y_k>=0)
    c = 0;
elseif(y_k<0 && y_k>=-h)
    c = c_max*(3*y_k^2/h^2+2*y_k^3/h^3);
elseif(y_k<-h)
    c = c_max;
end
if(y_k<=0)
    Ry = max([k*(-y_k)^p-c*dy_k,0]);
else
    Ry = 0;
end
Rx = Cuolumb_friction(Ry,dx_k,v_s,v_d,u_s,u_d);
R = [Rx;Ry];
end
