%% Funkcja zwracająca macierz rotacji względem osi y
function [Ry] = RotY(q)
    Ry = [cos(q) 0 sin(q);
           0 1 0;
           -sin(q) 0 cos(q)];
end

