%% Funkcja STEP z oprogramowania ADAMS
function s = STEP(x,x_1,h_1,x_2,h_2)
if(x <= x_1)
    s = h_1;
elseif (x>x_1 && x<=x_2)
    s = h_1+(h_2-h_1)*(x-x_1)^2/(x_2-x_1)^2*(3-2*(x-x_1)/(x_2-x_1));
else
    s = h_2;
end
end