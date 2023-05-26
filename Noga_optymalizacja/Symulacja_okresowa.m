%% Symulacja nogi dla periodycznych skoków
%% Dane
t_max = tmax;
w_max = wmax;
w_c = wc;
dalfa_max = 4/4*dalfa(1);
dbeta_min = 4/4*dbeta(1);
m1 = m(1);
m2 = m(2);
m3 = m(3);
l_1 = l(1);
l_2 = l(2);
M = Mc;
% Ilość skoków:
n = 10;
%% Wektory startowe:
[w,n_start] = min(abs(dalfa));
alfa_start = zeros(N1-n_start+1,1);
beta_start = zeros(N1-n_start+1,1);
dalfa_start = zeros(N1-n_start+1,1);
dbeta_start = zeros(N1-n_start+1,1);
t1_start = zeros(N1-n_start+1,1);
t2_start = zeros(N1-n_start+1,1);
alfa_start(1:N1-n_start+1) = alfa(n_start:N1);
dalfa_start(1:N1-n_start+1) = dalfa(n_start:N1);
beta_start(1:N1-n_start+1) = beta(n_start:N1);
dbeta_start(1:N1-n_start+1) = dbeta(n_start:N1);
t1_start(1:N1-n_start+1) = t1(n_start:N1);
t2_start(1:N1-n_start+1) = t2(n_start:N1);
time_start_N = linspace(0,time_N(N1)-time_N(n_start),N1-n_start+1)';
%% Wektory dla lotu:
[alfa_j,beta_j,dalfa_j,dbeta_j,t_j,N2] = jump_polynomials(m,l,N1,T,dy_,dy(N1),alfa_,dalfa_,alfa(N1),dalfa(N1),beta_,dbeta_,beta(N1),dbeta(N1));
alfa_lot = alfa_j;
dalfa_lot = dalfa_j;
beta_lot = beta_j;
dbeta_lot = dbeta_j;
%% Symulacja pierwszego wyskoku:

% Warunki początkowe:
q0 = [y(n_start);alfa(n_start);beta(n_start)];
dq0 = [dy(n_start);dalfa(n_start);dbeta(n_start)];
Y0 = [q0;dq0];
t0 = 0;
tK = time_N(N1)-time_N(n_start);
[time_start,Y_start] = ode45(@(time_start,Y_start) H_start(time_start,Y_start,time_start_N,alfa_start,dalfa_start,beta_start,dbeta_start,t1_start,t2_start,l,gamma,t_max,w_max,w_c), [t0,tK],Y0,odeset('RelTol',1e-6,'AbsTol',1e-6));
time_f = time_start;
Y = Y_start;
time_lot_K = t_j;
time_event = 0;
%% Symulacja dla kolejnych skoków:
for i = 1:n
    %% Symulacja dla lotu:
    time_lot_N = linspace(time_f(end),time_f(end)+time_lot_K,N2-1)';
    [time_lot,Y_lot] = ode45(@(time_lot,Y_lot) H_lot(time_lot,Y_lot,time_lot_N,alfa_lot,dalfa_lot,beta_lot,dbeta_lot,l,gamma,t_max,w_max,w_c), [time_f(end),time_f(end)+time_lot_K],Y(end,:),odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@(time_lot2,Y_lot2) collision(time_lot2,Y_lot2,dalfa_max,dbeta_min)));
    time_f = [time_f;time_lot];
    Y = [Y;Y_lot];
    alfa_lot2 = Y(end,2);
    beta_lot2 = Y(end,3);
    dalfa_lot2 = Y(end,5);
    dbeta_lot2 = Y(end,6);
     %% Symulacja dla lotu (dodatkowe oczekiwanie na upadek):
    [time_lot2,Y_lot2,te] = ode45(@(time_lot2,Y_lot2) H_lot2(time_lot2,Y_lot2,alfa_lot2,dalfa_lot2,beta_lot2,dbeta_lot2,l,gamma,t_max,w_max,w_c), [time_f(end),time_f(end)+2],Y(end,:),odeset('RelTol',1e-6,'AbsTol',1e-6,'Events',@(time_lot2,Y_lot2) collision(time_lot2,Y_lot2,dalfa_max,dbeta_min)));
    time_f = [time_f;time_lot2];
    Y = [Y;Y_lot2];
    time_event = [time_event;te];
    %% Symulacja dla upadku:
    time_upadek_N = linspace(time_f(end),time_f(end)+T,N1)';
    [time_upadek,Y_upadek] = ode45(@(time_upadek,Y_upadek)  H_upadek(time_upadek,Y_upadek,time_upadek_N,alfa,dalfa,beta,dbeta,t1,t2,l,gamma,t_max,w_max,w_c), [time_f(end),time_f(end)+T],Y(end,:),odeset('RelTol',1e-6,'AbsTol',1e-6));
    time_f = [time_f;time_upadek];
    Y = [Y;Y_upadek];
end
%% Wyniki:
% Wyniki symulacji
time_event = time_event(2:end);
t = time_f;
qs = Y(:,1:3);
dqs = Y(:,4:6);
ys = qs(:,1);
alfas = qs(:,2);
betas = qs(:,3);
dys = dqs(:,1);
dalfas = dqs(:,2);
dbetas = dqs(:,3);

% Generowanie momentów napędowych, sił reakcji podłoża oraz trajekotrii referencyjnych
[alfa_s,dalfa_s,beta_s,dbeta_s,T_s,time_S,N2] = gen_trajectory(m,l,N1,T,alfa,dalfa,beta,dbeta,dy,dy_,alfa_,dalfa_,beta_,dbeta_,n);
n_c  = size(t);
Rs = zeros(n_c(1),2); % Siły reakcji i tarcia podłoża z symulacji
n_e = size(time_event);
t1s = zeros(n_c); % Moment napędowy 1 z symulacji 
t2s = zeros(n_c); % Moment napędowy 2 z symulacji 
Fy_op = zeros(n_c); % Siły reakcji podłoża referencyjne
Fx_op = zeros(n_c); % Siły tarcia podłoża referencyjne
n_s = size(time_start);
t1s(1:n_s(1),1) = interp1(time_N(1:N1-n_start+1),t1(n_start:N1),time_start);
t2s(1:n_s(1),1) = interp1(time_N(1:N1-n_start+1),t2(n_start:N1),time_start);
Fy_op(1:n_s(1),1) = interp1(time_N(1:N1-n_start+1),Fy(n_start:N1),time_start);
Fx_op(1:n_s(1),1) = interp1(time_N(1:N1-n_start+1),Fx(n_start:N1),time_start);
for i = 1:n_e(1)
    [w,N_e1] = min(abs(time_f-time_event(i)));
    [w,N_e2] = min(abs(time_f-time_event(i)-T));
    time_O = linspace(time_event(i),time_event(i)+T,N1)';
    
    t1s(N_e1:N_e2) = interp1(time_O,t1,time_f(N_e1:N_e2));
    t2s(N_e1:N_e2) = interp1(time_O,t2,time_f(N_e1:N_e2));
    Fy_op(N_e1:N_e2) = interp1(time_O,Fy,time_f(N_e1:N_e2));
    Fx_op(N_e1:N_e2) = interp1(time_O,Fx,time_f(N_e1:N_e2));
end
% Generowanie współrzędnych i prędkości złączowych referencyjnych
alfa_r = interp1(time_S,alfa_s,time_f); 
dalfa_r = interp1(time_S,dalfa_s,time_f);
beta_r = interp1(time_S,beta_s,time_f);
dbeta_r = interp1(time_S,dbeta_s,time_f);
y_com = zeros(n_c(1),1); % Położenie środka masy z symulacji
dy_com = zeros(n_c(1),1); % Prędkość środka masy z symulacji
for i = 1:n_c(1)
    Rs(i,:) = reaction_force_spring([ys(i);alfas(i);betas(i)],[dys(i);dalfas(i);dbetas(i)],l)';
    tn = PD_cart_regulator([alfa_r(i);beta_r(i)],[alfas(i);betas(i)],[dalfa_r(i);dbeta_r(i)],[dalfas(i);dbetas(i)],l);
    tn = tn + PD_joint_regulator([alfa_r(i);beta_r(i)],[alfas(i);betas(i)],[dalfa_r(i);dbeta_r(i)],[dalfas(i);dbetas(i)]);
    t1s(i) = t1s(i) + tn(1);
    t2s(i) = t2s(i) + tn(2);    
    if(t1s(i) > gamma(1)*t_max(1))
      t1s(i) = gamma(1)*t_max(1);
    elseif(t1s(i)<-gamma(1)*t_max(1))
      t1s(i) = -gamma(1)*t_max(1);
    end
    if(t2s(i) > gamma(2)*t_max(2))
      t2s(i) = gamma(2)*t_max(2);
    elseif(t2s(i)<-gamma(2)*t_max(2))
     t2s(i) = -gamma(2)*t_max(2);
    end
    t1s(i) = BLDC_motor_constraint(t1s(i),dalfas(i),gamma(1),t_max(1),w_max(1),w_c(1));
    t2s(i) = BLDC_motor_constraint(t2s(i),dbetas(i),gamma(2),t_max(2),w_max(2),w_c(2));
    y_com(i) = (ys(i)*m1+(ys(i)+l_1/2*sin(alfas(i)))*m2+(ys(i)+l_1*sin(alfas(i))+l_2/2*sin(alfas(i)+betas(i))*m3))/(M);
    dy_com(i) = (dys(i)*m1+(dys(i)+l_1/2*cos(alfas(i))*dalfas(i))*m2+(dys(i)+l_1*cos(alfas(i))*dalfas(i)+l_2/2*cos(alfas(i)+betas(i))*(dalfas(i)+dbetas(i))*m3))/(M);
end

%% Funkcja dla startu:
function dY_start = H_start(time_start,Y_start,time_start_N,alfa_start,dalfa_start,beta_start,dbeta_start,t1_start,t2_start,l,gamma,t_max,w_max,w_c)
q = Y_start(1:3,:);
dq = Y_start(4:6,:);
tn1 = interp1(time_start_N,t1_start,time_start);
tn2 = interp1(time_start_N,t2_start,time_start);
alfa_r = interp1(time_start_N,alfa_start,time_start);
dalfa_r = interp1(time_start_N,dalfa_start,time_start);
beta_r = interp1(time_start_N,beta_start,time_start);
dbeta_r = interp1(time_start_N,dbeta_start,time_start);
tn = [tn1;tn2];
tn = tn + PD_cart_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1),l);
tn = tn + PD_joint_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1));
if(tn(1) > gamma(1)*t_max(1))
   tn(1) = gamma(1)*t_max(1);
elseif(tn(1)<-gamma(1)*t_max(1))
   tn(1) = -gamma(1)*t_max(1);
end
if(tn(2) > gamma(2)*t_max(2))
   tn(2) = gamma(2)*t_max(2);
elseif(tn(2)<-gamma(2)*t_max(2))
  tn(2) = -gamma(2)*t_max(2);
end
tn(1) = BLDC_motor_constraint(tn(1),dq(2),gamma(1),t_max(1),w_max(1),w_c(1));
tn(2) = BLDC_motor_constraint(tn(2),dq(3),gamma(2),t_max(2),w_max(2),w_c(2));
tnc = [0;tn(1);tn(2)];
dY_start = dyn_simulation_noga_2D_spring(q,dq,tnc,l);
end

%% Funkcja dla lotu:
function dY_lot = H_lot(time_lot,Y_lot,time_lot_N,alfa_lot,dalfa_lot,beta_lot,dbeta_lot,l,gamma,t_max,w_max,w_c)
q = Y_lot(1:3,:);
dq = Y_lot(4:6,:);
alfa_r = interp1(time_lot_N,alfa_lot,time_lot);
dalfa_r = interp1(time_lot_N,dalfa_lot,time_lot);
beta_r = interp1(time_lot_N,beta_lot,time_lot);
dbeta_r = interp1(time_lot_N,dbeta_lot,time_lot);
tn = [0;0];
tn = tn + PD_cart_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1),l);
tn = tn + PD_joint_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1));
if(tn(1) > gamma(1)*t_max(1))
   tn(1) = gamma(1)*t_max(1);
elseif(tn(1)<-gamma(1)*t_max(1))
   tn(1) = -gamma(1)*t_max(1);
end
if(tn(2) > gamma(2)*t_max(2))
   tn(2) = gamma(2)*t_max(2);
elseif(tn(2)<-gamma(2)*t_max(2))
  tn(2) = -gamma(2)*t_max(2);
end
tn(1) = BLDC_motor_constraint(tn(1),dq(2),gamma(1),t_max(1),w_max(1),w_c(1));
tn(2) = BLDC_motor_constraint(tn(2),dq(3),gamma(2),t_max(2),w_max(2),w_c(2));
tnc = [0;tn(1);tn(2)];
dY_lot = dyn_simulation_noga_2D_spring(q,dq,tnc,l);
end

%% Funkcja dla lotu 2 (gdy upadek się wydłuża):
function dY_lot2 = H_lot2(time_lot2,Y_lot2,alfa_lot2,dalfa_lot2,beta_lot2,dbeta_lot2,l,gamma,t_max,w_max,w_c)
q = Y_lot2(1:3,:);
dq = Y_lot2(4:6,:);
alfa_r = alfa_lot2;
dalfa_r = dalfa_lot2;
beta_r = beta_lot2;
dbeta_r = dbeta_lot2;
tn = [0;0];
tn = tn + PD_cart_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1),l);
tn = tn + PD_joint_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1));
if(tn(1) > gamma(1)*t_max(1))
   tn(1) = gamma(1)*t_max(1);
elseif(tn(1)<-gamma(1)*t_max(1))
   tn(1) = -gamma(1)*t_max(1);
end
if(tn(2) > gamma(2)*t_max(2))
   tn(2) = gamma(2)*t_max(2);
elseif(tn(2)<-gamma(2)*t_max(2))
  tn(2) = -gamma(2)*t_max(2);
end
tn(1) = BLDC_motor_constraint(tn(1),dq(2),gamma(1),t_max(1),w_max(1),w_c(1));
tn(2) = BLDC_motor_constraint(tn(2),dq(3),gamma(2),t_max(2),w_max(2),w_c(2));
tnc = [0;tn(1);tn(2)];
dY_lot2 = dyn_simulation_noga_2D_spring(q,dq,tnc,l);
end

%% Funkcja dla upadku:
function dY_upadek = H_upadek(time_upadek,Y_upadek,time_upadek_N,alfa,dalfa,beta,dbeta,t1,t2,l,gamma,t_max,w_max,w_c)
q = Y_upadek(1:3,:);
dq = Y_upadek(4:6,:);
tn1 = interp1(time_upadek_N,t1,time_upadek);
tn2 = interp1(time_upadek_N,t2,time_upadek);
alfa_r = interp1(time_upadek_N,alfa,time_upadek);
dalfa_r = interp1(time_upadek_N,dalfa,time_upadek);
beta_r = interp1(time_upadek_N,beta,time_upadek);
dbeta_r = interp1(time_upadek_N,dbeta,time_upadek);
tn = [tn1;tn2];
tn = tn + PD_cart_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1),l);
tn = tn + PD_joint_regulator([alfa_r;beta_r],q(2:3,1),[dalfa_r;dbeta_r],dq(2:3,1));
if(tn(1) > gamma(1)*t_max(1))
   tn(1) = gamma(1)*t_max(1);
elseif(tn(1)<-gamma(1)*t_max(1))
   tn(1) = -gamma(1)*t_max(1);
end
if(tn(2) > gamma(2)*t_max(2))
   tn(2) = gamma(2)*t_max(2);
elseif(tn(2)<-gamma(2)*t_max(2))
  tn(2) = -gamma(2)*t_max(2);
end
tn(1) = BLDC_motor_constraint(tn(1),dq(2),gamma(1),t_max(1),w_max(1),w_c(1));
tn(2) = BLDC_motor_constraint(tn(2),dq(3),gamma(2),t_max(2),w_max(2),w_c(2));
tnc = [0;tn(1);tn(2)];
dY_upadek = dyn_simulation_noga_2D_spring(q,dq,tnc,l);
end

%% Funkcja event (sprawdza spełnienie warunków upadku):
function [value,isterminal,direction] = collision(time_lot,Y_lot,dalfa_max,dbeta_min)
dq = Y_lot(4:6,:);
value = [dq(2)-dalfa_max; dq(3)-dbeta_min];  
isterminal = [1;  1];         % stop at local minimum
direction  = [1; -1];         % [local minimum, local maximum]
end