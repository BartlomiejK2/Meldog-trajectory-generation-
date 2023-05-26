%% Nieliniowe równości i nierówności optymalizacji
function [c,ceq] = nonlcon_noga_2D(x,N,l,y0,M,m)
%% Stałe:
g = 9.81;
m1 = m(1);
m2 = m(2);
m3 = m(3);
l_1 = l(1);
l_2 = l(2);
%% Inicjalizacja zmiennych:
y_ = x(1);
y = x(2:N+1);
alfa_ = x(N+2);
alfa = x(N+3:2*N+2);
beta_ = x(2*N+3);
beta = x(2*N+4:3*N+3);
dy_ = x(3*N+4);
dy = x(3*N+5:4*N+4);
dalfa_ = x(4*N+5);
dalfa = x(4*N+6:5*N+5);
dbeta_ = x(5*N+6);
dbeta = x(5*N+7:6*N+6);
t1 = x(6*N+7:7*N+6);
t2 = x(7*N+7:8*N+6);
Fx = x(8*N+7:9*N+6);
Fy = x(9*N+7:10*N+6);
T = x(10*N+7);
p = x(10*N+8);
Rdt = x(10*N+9:10*N+10);
dt = T/(N-1);
f1 = zeros(6*(N-1),1);
f2 = zeros(2*N,1);
f3 = zeros(2*N,1);

%% Dynamika:
for i = 1:N-1
    q_i = [y(i);alfa(i);beta(i)];
    q_i_1 = [y(i+1);alfa(i+1);beta(i+1)];
    dq_i = [dy(i);dalfa(i);dbeta(i)];
    dq_i_1 = [dy(i+1);dalfa(i+1);dbeta(i+1)];
    t_i = [0;t1(i);t2(i)];
    t_i_1 = [0;t1(i+1);t2(i+1)];
    F_i = [Fx(i);Fy(i);0];
    F_i_1 = [Fx(i+1);Fy(i+1);0];
    f1(6*(i-1)+1:6*i,1) = [q_i_1;dq_i_1] - [q_i;dq_i] - dt/2*(forward_dyn_noga_2D(q_i,dq_i,t_i,F_i,zeros(3,1)) + forward_dyn_noga_2D(q_i_1,dq_i_1,t_i_1,F_i_1,zeros(3,1)));
end

%% Kinematyka:
for j = 1:N
    q_j = [y(j);alfa(j);beta(j)];
    dq_j = [dy(j);dalfa(j);dbeta(j)];
    f2(2*(j-1)+1:2*j,1) = forward_kin_noga_2D(q_j,l) - [p;0];
    f3(2*(j-1)+1:2*j,1) =  forward_vel_noga_2D(q_j,dq_j,l);
end

%% Warunki upadku:
q0 = [y(1);alfa(1);beta(1)];
dq0 = [dy(1);dalfa(1);dbeta(1)];
dq_ = [dy_;dalfa_;dbeta_];

%y i dy N:
y_n = y(N);
alfa_n = alfa(N);
beta_n = beta(N);
dy_n = dy(N);
dalfa_n = dalfa(N);
dbeta_n = dbeta(N);
y_c = (y_n*m1+(y_n+l_1/2*sin(alfa_n))*m2+(y_n+l_1*sin(alfa_n)+l_2/2*sin(alfa_n+beta_n)*m3))/(M);
dy_c = (dy_n*m1+(dy_n+l_1/2*cos(alfa_n)*dalfa_n)*m2+(dy_n+l_1*cos(alfa_n)*dalfa_n+l_2/2*cos(alfa_n+beta_n)*(dalfa_n+dbeta_n)*m3))/(M);

%y i dy _ :
y_c_ = (y_*m1+(y_+l_1/2*sin(alfa_))*m2+(y_+l_1*sin(alfa_)+l_2/2*sin(alfa_+beta_)*m3))/(M);
dy_c_ = (dy_*m1+(dy_+l_1/2*cos(alfa_)*dalfa_)*m2+(dy_+l_1*cos(alfa_)*dalfa_+l_2/2*cos(alfa_+beta_)*(dalfa_+dbeta_)*m3))/(M);

f4 = [(dy_c_)^2-(dy_c)^2 - 2*g*(y_c-y_c_);
    plastic_impact(q0,dq0,dq_,Rdt,l)];

    
    
ceq = [f1;f2;f3;f4];
%% Nierówności:

%% Prędkość środka masy:
dy_c_max = 2.5;
c = -dy_c + dy_c_max; 
end

