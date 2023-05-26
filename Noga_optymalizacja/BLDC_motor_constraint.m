%% Funkcja generująca moment napędowy silnika BLDC zgodnie z jego charakterystyką
function t = BLDC_motor_constraint(t_r,w_r,g,t_max,w_max,w_c)
b = (g*t_max)/(w_max-w_c);
A = [1 g*b;
    -1 g*b;
    1 -g*b;
    -1 -g*b];
b_r = [b*w_max;
    b*w_max;
    b*w_max;
    b*w_max];
x = A*[t_r;
    w_r]-b_r;
if(x(1)>0)
    t = b*w_max-g*b*w_r;
elseif(x(2)>0)
    t = -(b*w_max-g*b*w_r);
elseif(x(3)>0)
    t = (b*w_max+g*b*w_r);
elseif(x(4)>0)
    t = -(b*w_max+g*b*w_r);
else
    t = t_r;
end
if(w_r>w_max || w_r<-w_max)
    t = 0;
end
end

