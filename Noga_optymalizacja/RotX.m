%% Funkcja zwracająca macierz rotacji względem osi x
function [Rx] = RotX(q)
    Rx = [1 0 0;
           0 cos(q) -sin(q);
           0 sin(q) cos(q)];
end

