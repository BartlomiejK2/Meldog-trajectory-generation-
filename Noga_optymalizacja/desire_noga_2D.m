%% Funkcja celu optymalizacji
function [f] = desire_noga_2D(x,N,l)
alfa = 1;
beta = 0.4;
l_1 = l(1);
l_2 = l(2);
m1 = 2;
m2 = 0.5;
m3 = 0.5;
y_n = x(N+1,1);
alfa_n = x(2*N+2,1);
beta_n = x(3*N+3,1);
dy_n = x(4*N+4,1);
dalfa_n = x(5*N+5,1);
dbeta_n = x(6*N+6,1);
y_c = (y_n*m1+(y_n+l_1/2*sin(alfa_n))*m2+(y_n+l_1*sin(alfa_n)+l_2/2*sin(alfa_n+beta_n)*m3))/(m1+m2+m3);
dy_c = (dy_n*m1+(dy_n+l_1/2*cos(alfa_n)*dalfa_n)*m2+(dy_n+l_1*cos(alfa_n)*dalfa_n+l_2/2*cos(alfa_n+beta_n)*(dalfa_n+dbeta_n)*m3))/(m1+m2+m3);
g = 9.81;
h = y_c +1/2*(dy_c)^2/g; % Maksymalna wysokość środka masy
Rdt = x(10*N+9:10*N+10);
n_Rdt = Rdt'*Rdt; % Kwadrat normy impulsu
T = x(10*N+7);
dt = T/(N-1);
t1 = x(6*N+7:7*N+6);
t2 = x(7*N+7:8*N+6);
St = 0;
for i = 1:N-1
    St = St+1/2*(t1(i+1)^2 + t2(i+1)^2 + t1(i)^2 + t2(i)^2)*dt; % Całka kwadratu momentu napędowego
end
f = St;
end

