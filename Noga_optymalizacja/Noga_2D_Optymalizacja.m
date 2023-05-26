%% SKRYPT DO OPTYMALIZACJI RUCHU NOGI 2D

%% Budowa wektora x:
%{
y_ = x(1); - Wysokość członu 1 chwile przed upadkiem
y = x(2:N1+1); - Wysokość członu 1
alfa_ = x(N1+2); - Kąt obrotu w pierwszej parze obrotowej chwile przed upadkiem
alfa = x(N1+3:2*N1+2); - Kąt obrotu w pierwszej parze obrotowej
beta_ = x(2*N1+3); - Kąt obrotu w drugiej parze obrotowej chwile przed upadkiem
beta = x(2*N1+4:3*N1+3); - Kąt obrotu w drugiej parze obrotowej 
dy_ = x(3*N1+4); - Prędkość liniowa członu 1 chwile przed upadkiem
dy = x(3*N1+5:4*N1+4); - Prędkość liniowa członu 1 
dalfa_ = x(4*N1+5); - Prędkość kątowa w pierwszej parze obrotowej chwile przed upadkiem
dalfa = x(4*N1+6:5*N1+5); - Prędkość kątowa w pierwszej parze obrotowej
dbeta_ = x(5*N1+6); - Prędkość kątowa w drugiej parze obrotowej chwile przed upadkiem
dbeta = x(5*N1+7:6*N1+6); - Prędkość kątowa w drugiej parze obrotowej 
t1 = x(6*N1+7:7*N1+6); - Moment napędowy w pierwszej parze obrotowej
t2 = x(7*N1+7:8*N1+6); - Moment napędowy w drugiej parze obrotowej
Fx = x(8*N1+7:9*N1+6); - Siła reakcji podłoża
Fy = x(9*N1+7:10*N1+6); - Siła tarcia podłoża
T = x(10*N1+7); - Czas kontaktu z podłożem
p = x(10*N1+8); - Odległość końcówki nogi względem układu 0
Rdt = x(10*N1+9:10*N1+10); - Wartość impaktu
%}
%% Opcje
clear all
options = optimoptions('fmincon','MaxFunEvals',700000,'MaxIter',2500,'PlotFcn','optimplotconstrviolation','StepTolerance',1e-10,'Display','iter');

%% Dane modelu
u = 0.7;
tmax = [2.5;2.5];
gamma = [16;16];
wmax = [60*2*3.14;60*2*3.14];
wc = [45*2*3.14;45*2*3.14];
l = [0.25;0.25];
Mc = 3;
m = [2;0.5;0.5];
%% 1 Optymalizacja 
N1 = 30;

% Nierówności
[ub,lb] = ub_lb_noga_2D(N1,tmax,wmax,l,gamma);
% Równości równania:
[Aeq,beq] =  eq_noga_2D(N1);
% Nierówności równiania:
[A,b] = neq_noga_2D(N1,u,gamma,tmax,wc,wmax);
% Warunki początkowe:
y0 = zeros(N1+1,1);
alfa0 = zeros(N1+1,1);
beta0 = zeros(N1+1,1);
dy0 = zeros(N1+1,1);
dalfa0 = zeros(N1+1,1);
dbeta0 = zeros(N1+1,1);
t10 = zeros(N1,1);
t20 = zeros(N1,1);
Fx0 = zeros(N1,1);
Fy0 = zeros(N1,1);
T0 = 0.3;
p0 = 0;
y_0 = 0.1;
% Upadek i skok:
y_min = 0.1;
y_max = 0.40;
b_us = [y_0;y_min;y_max];
A_us = [1 1 1;
    N1^2/9 N1/3 1;
    N1^2 N1 1];
lambda_us = A_us\b_us;
% Sam wyskok:
y_N = 0.4;
ay = (y_N-y_0)/(N1^2-2*N1+1);
by = -2*ay;
cy = ay+y_0;
%% Generowanie współrzędnych złączowych startowych:
for i = 2:N1+1
    %y0(i) = ay*i^2+by*i+cy; % Sam wyskok
    y0(i) = lambda_us (1)*i^2+lambda_us (2)*i+lambda_us (3); % Upadek i skok
    % Odwrotna kinematyka położenia dla warunków początkowych
    inv_h = inv_kin_noga_2D(y0(i),p0,l);
    alfa0(i) =  inv_h(1);
    beta0(i) =  inv_h(2);
end
%% Generowanie prędkości złączowych startowych:
for i = 2:N1+1
    %dy0(i) = (2*ay*i+by)*N1/T0; % Sam wyskok
    dy0(i) = (2*lambda_us (1)*i+lambda_us (2))*N1/T0; % Upadek i skok
    % Odwrotna kinematyka prędkości dla warunków początkowych
    inv_dh = inv_vel_noga_2D([y0(i);alfa0(i);beta0(i)],dy0(i),l);
    dalfa0(i) = inv_dh(1);
    dbeta0(i) = inv_dh(2);
end
%% Generowanie momentów napędowych oraz siły reakcji podłoża startowych:
for i = 1:N1
    Fy0(i) = 150;
    t10(i) = -tmax(1)*gamma(1)/2;
    t20(i) = tmax(2)*gamma(2)/2;
end
y0(1) = y0(2);
alfa0(1) = alfa0(2);
beta0(1) = beta0(2);
dy0(1) = dy0(2);
dalfa0(1) = dalfa0(2);
dbeta0(1) = dbeta0(2);
Fy0(N1) = 0;
q0 = [y0;alfa0;beta0];
dq0 = [dy0;dalfa0;dbeta0];
t0 = [t10;t20];
F0 = [Fx0;Fy0];
Rdt0 = zeros(2,1);
x0 = [q0;dq0;t0;F0;T0;p0;Rdt0];

%% Solver:
nonlcon = @(x) nonlcon_noga_2D(x,N1,l,y_0,Mc,m);
fun = @(x) desire_noga_2D(x,N1,l);
tic
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);
toc

%% Post-solving:
y_ = x(1);
y = x(2:N1+1);
alfa_ = x(N1+2);
alfa = x(N1+3:2*N1+2);
beta_ = x(2*N1+3);
beta = x(2*N1+4:3*N1+3);
dy_ = x(3*N1+4);
dy = x(3*N1+5:4*N1+4);
dalfa_ = x(4*N1+5);
dalfa = x(4*N1+6:5*N1+5);
dbeta_ = x(5*N1+6);
dbeta = x(5*N1+7:6*N1+6);
t1 = x(6*N1+7:7*N1+6);
t2 = x(7*N1+7:8*N1+6);
Fx = x(8*N1+7:9*N1+6);
Fy = x(9*N1+7:10*N1+6);
T = x(10*N1+7);
p = x(10*N1+8);
Rdt = x(10*N1+9:10*N1+10);
time_N = linspace(0,T,N1);
f = desire_noga_2D(x,N1,l);
[c,ceq] = nonlcon_noga_2D(x,N1,l,y_0,Mc,m);
n_ceq = norm(ceq);
n_eq = norm(Aeq*x-beq);
f0 = desire_noga_2D(x0,N1,l);
[c0,ceq0] = nonlcon_noga_2D(x0,N1,l,y_0,Mc,m);
n_ceq0 = norm(ceq0);
n_eq0 = norm(Aeq*x0-beq);
h = calc_height(x,N1,l);
