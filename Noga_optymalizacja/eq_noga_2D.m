%% Równości optymalizacji
function [Aeq,beq] = eq_noga_2D(N)

Aeq = zeros(9,10*N+10);
beq = zeros(9,1);
%% Warunki wyskoku:

% Wektor sił reakcji podłoża:
Aeq(1,9*N+6) = 1;
Aeq(2,10*N+6) = 1;
%% Warunki upadku:

% Położenie:
Aeq(3,0*(N+1)+2) = 1;
Aeq(3,0*(N+1)+1) = -1;

Aeq(4,1*(N+1)+2) = 1;
Aeq(4,1*(N+1)+1) = -1;

Aeq(5,2*(N+1)+2) = 1;
Aeq(5,2*(N+1)+1) = -1;

% Prędkości:
Aeq(6,4*(N+1)+1) = 1;
Aeq(7,5*(N+1)+1) = 1;
end



