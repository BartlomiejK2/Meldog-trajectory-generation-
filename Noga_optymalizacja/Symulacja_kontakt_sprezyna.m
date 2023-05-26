%% Skrypt rozwiązujący symulacje nogi w trakcie kontaku z podłożem, zamodelowanym jako nieliniowa sprężyna

%% Przyotowanie symulacji:
N1l = linspace(0,T,N1)';
q0 = [y(1);alfa(1);beta(1)];
dq0 = [dy(1);dalfa(1);dbeta(1)];
Y0 = [q0;dq0];
t0 = 0;
tK = T;
%% Symulacja:
[t,Y] = ode45(@(t,Y) H(t,Y,N1l,t1,t2,l), [t0,tK],Y0,odeset('RelTol',1e-8,'AbsTol',1e-8));
%% Wyniki:
q = Y(:,1:3); % Współrzędne złączowe
dq = Y(:,4:6); % Prędkości złączowe
ys = q(:,1);
alfas = q(:,2);
betas = q(:,3);
dys = dq(:,1);
dalfas = dq(:,2);
dbetas = dq(:,3);
t1s = interp1(N1l,t1,t); % Moment napędowy 1
t2s = interp1(N1l,t2,t); % Moment napędowy 2
n = size(t);
Rs = zeros(n(1),2); % Siła reakcji i tarcia podłoża 
for i = 1:n(1)
    Rs(i,:) = reaction_force_spring([ys(i);alfas(i);betas(i)],[dys(i);dalfas(i);dbetas(i)],l)';
end
%% Równania stanu:
function dY = H(t,Y,N1l,t1,t2,l)
t1 = interp1(N1l,t1,t);
t2 = interp1(N1l,t2,t);
t = [0;t1;t2];
q = Y(1:3,:);
dq = Y(4:6,:);
dY = dyn_simulation_noga_2D_spring(q,dq,t,l);
end