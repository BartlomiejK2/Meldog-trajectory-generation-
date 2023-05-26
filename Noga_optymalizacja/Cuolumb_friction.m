%% Funkcja generująca siłę tarcia według metody Cuolumba zawartej w oprogramowaniu ADAMS
function R_f = Cuolumb_friction(Fy,v,v_s,v_d,u_s,u_d)
if(v == -v_s)
    u = u_s;
elseif(v == v_s)
    u = -u_s;
elseif(v == 0)
    u = 0;
elseif(v == -v_d)
    u = u_d;
elseif(v == v_d)
    u = -u_d;
elseif(abs(v)>v_d)
    u = -sign(v) * u_d;
elseif(abs(v)>= v_s && abs(v)<=v_d)
    u = -STEP(abs(v),v_s,u_s,v_d,u_d)*sign(v);
elseif(v>= -v_s && v<=v_s)
    u = STEP(v,-v_s,u_s,v_s,-u_s);
end
R_f = u*Fy;
end
function s = STEP(x,x_1,h_1,x_2,h_2)
if(x <= x_1)
    s = h_1;
elseif (x>x_1 && x<=x_2)
    s = h_1+(h_2-h_1)*(x-x_1)^2/(x_2-x_1)^2*(3-2*(x-x_1)/(x_2-x_1));
else
    s = h_2;
end
end