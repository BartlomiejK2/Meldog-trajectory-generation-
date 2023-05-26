%% Ograniczenia optymalizacji
function [ub,lb] = ub_lb_noga_2D(N,tmax,wmax,l,g)
l_sum = l(1)+l(2);
ub = zeros(10*N+10,1);
lb = zeros(10*N+10,1);
% y
ub(1:N+1) = l_sum;
lb(1:N+1) = 0.05;

% alfa
ub((N+1)+1:2*(N+1),1) = pi/2;
lb((N+1)+1:2*(N+1),1) = -pi/2;

% beta
ub(2*(N+1)+1:3*(N+1),1) = 0;
lb(2*(N+1)+1:3*(N+1),1) = -5/6*pi;

% dy
ub(3*(N+1)+1:4*(N+1),1) = Inf;
lb(3*(N+1)+1:4*(N+1),1) = -Inf;

% dalfa
ub(4*(N+1)+1:5*(N+1),1) = wmax(1)/g(1);
lb(4*(N+1)+1:5*(N+1),1) = -wmax(1)/g(1);

% dbeta
ub(5*(N+1)+1:6*(N+1),1) = wmax(2)/g(2);
lb(5*(N+1)+1:6*(N+1),1) = -wmax(2)/g(2);

% t1
ub(6*(N+1)+1:7*N+6,1) = g(1)*tmax(1);
lb(6*(N+1)+1:7*N+6,1) = -g(1)*tmax(1);

% t2
ub(7*N+6+1:8*N+6,1) = g(2)*tmax(2);
lb(7*N+6+1:8*N+6,1) = -g(2)*tmax(2);

% Fx
ub(8*N+6+1:9*N+6,1) = Inf;
lb(8*N+6+1:9*N+6,1) = -Inf;

% Fy
ub(9*N+6+1:10*N+6,1) = 250;
lb(9*N+6+1:10*N+6,1) = 0;

% T
ub(10*N+6+1,1) = Inf;
lb(10*N+6+1,1) = 0.2;

% p
ub(10*N+6+2,1) = l_sum/4;
lb(10*N+6+2,1) = -l_sum/4;

%Rdt
ub(10*N+6+3,1) = Inf;
lb(10*N+6+3,1) = -Inf;
ub(10*N+6+4,1) = Inf;
lb(10*N+6+4,1) = 0;

% Warunek wyskoku
%lb(4*(N+1),1) = 3;
%lb(1*N,1) = 0.40*l_sum;

end


