%% Funkcja do generowania trajektorii referenfycjnych nogi 
function [alfa_s,dalfa_s,beta_s,dbeta_s,t_s,time_S,N2] = gen_trajectory(m,l,N1,T,alfa,dalfa,beta,dbeta,dy,dy_,alfa_,dalfa_,beta_,dbeta_,n)


%% Obliczanie punktu startowego oraz generacja alfy i bety w czasie skoku:
[w,n_start] = min(abs(dalfa));
[alfa_j,beta_j,dalfa_j,dbeta_j,t_j,N2] = jump_polynomials(m,l,N1,T,dy_,dy(N1),alfa_,dalfa_,alfa(N1),dalfa(N1),beta_,dbeta_,beta(N1),dbeta(N1));
%% Tworzenie wektorów:
alfa_s = zeros(N1-n_start+1+n*(N1+N2-1),1);
beta_s = zeros(N1-n_start+1+n*(N1+N2-1),1);
dalfa_s = zeros(N1-n_start+1+n*(N1+N2-1),1);
dbeta_s = zeros(N1-n_start+1+n*(N1+N2-1),1);
%% Generowanie całości
alfa_s(1:N1-n_start+1) = alfa(n_start:N1);
dalfa_s(1:N1-n_start+1) = dalfa(n_start:N1);
beta_s(1:N1-n_start+1) = beta(n_start:N1);
dbeta_s(1:N1-n_start+1) = dbeta(n_start:N1);
t_s = T*(N1-n_start+1)/N1;
for i = 1:n
    alfa_s(N1-n_start+2+(i-1)*(N1+N2-1):N1-n_start+2+(i)*(N1+N2-1)-1) = [alfa_j;alfa];
    dalfa_s(N1-n_start+2+(i-1)*(N1+N2-1):N1-n_start+2+(i)*(N1+N2-1)-1) = [dalfa_j;dalfa];
    beta_s(N1-n_start+2+(i-1)*(N1+N2-1):N1-n_start+2+(i)*(N1+N2-1)-1) = [beta_j;beta];
    dbeta_s(N1-n_start+2+(i-1)*(N1+N2-1):N1-n_start+2+(i)*(N1+N2-1)-1) = [dbeta_j;dbeta];
    t_s  = t_s + t_j + T;
end
time_S = linspace(0,t_s,N1-n_start+1+n*(N1+N2-1))';
end

