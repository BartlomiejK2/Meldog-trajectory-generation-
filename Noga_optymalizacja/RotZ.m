%% Funkcja zwracająca macierz rotacji względem osi z
function [Rz] = RotZ(q)
    Rz = [cos(q) -sin(q) 0;
           sin(q) cos(q) 0;
           0 0 1];
end

