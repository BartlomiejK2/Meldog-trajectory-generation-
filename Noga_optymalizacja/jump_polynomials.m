%% Funkcja generująca trajektorie dla lotu nogi
function [alfa_j,beta_j,dalfa_j,dbeta_j,t_j,N2] = jump_polynomials(m,l,N1,T,dy_0,dy_N,alfa_0,dalfa_0,alfa_N,dalfa_N,beta_0,dbeta_0,beta_N,dbeta_N)
m1 = m(1);
m2 = m(2);
m3 = m(3);
M = m1+m2+m3;
g = 9.81;
l_1 = l(1);
l_2 = l(2);
dyc_0 = (dy_0*m1+(dy_0+l_1/2*cos(alfa_0)*dalfa_0)*m2+(dy_0+l_1*cos(alfa_0)*dalfa_0+l_2/2*cos(alfa_0+beta_0)*(dalfa_0+dbeta_0)*m3))/(M);
dyc_N = (dy_N*m1+(dy_N+l_1/2*cos(alfa_N)*dalfa_N)*m2+(dy_N+l_1*cos(alfa_N)*dalfa_N+l_2/2*cos(alfa_N+beta_N)*(dalfa_N+dbeta_N)*m3))/(M);
tc = (abs(dyc_0)+abs(dyc_N))/g;
t_j = tc;
N2 = round(tc*N1/T)+1;
A = [0, 0, 0, 1;
    0, 0, 1, 0;
    tc^3, tc^2, tc, 1;
    3*tc^2, 2*tc, 1, 0];
b1 = [alfa_N,dalfa_N,alfa_0,dalfa_0]';
b2 = [beta_N,dbeta_N,beta_0,dbeta_0]';
x1 = A\b1;
x2 = A\b2;
dx1 = [3*x1(1);2*x1(2);x1(3)];
dx2 = [3*x2(1);2*x2(2);x2(3)];
time_N2 = linspace(0,tc,N2);
alfa_j = polyval(x1,time_N2)';
dalfa_j = polyval(dx1,time_N2)';
beta_j = polyval(x2,time_N2)';
dbeta_j = polyval(dx2,time_N2)';
alfa_j = alfa_j(2:N2,1);
dalfa_j = dalfa_j(2:N2,1);
beta_j = beta_j(2:N2,1);
dbeta_j = dbeta_j(2:N2,1);
end

