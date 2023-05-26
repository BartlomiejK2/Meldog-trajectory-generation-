%% Nierówności optymalizacji
function [A,b] = neq_noga_2D(N,u,g,tmax,wc,wmax)
A_F = zeros(3*N,10*N+10);
b_F = zeros(3*N,1);
A_t1 = zeros(4*N,10*N+10);
A_t2 = zeros(4*N,10*N+10);
b_t1 = zeros(4*N,1);
b_t2 = zeros(4*N,1);

%% Warunek sił reakcji podłoża (|Fx|<|Fy|)
A_s = [0, -1;
    1, -u;
    -1, -u];
for i = 1:N
    A_F(3*(i-1)+1:3*(i-1)+3,8*N+6+i) = A_s(1:3,1);
    A_F(3*(i-1)+1:3*(i-1)+3,9*N+6+i) = A_s(1:3,2);
end
%% Warunek zależnosci momentu od prędkości silnika BLDC: (|t|<=|g*b*|domega|+b|)
b1 = -g(1)*tmax(1)/(wc(1)-wmax(1));
b2 = -g(2)*tmax(2)/(wc(2)-wmax(2));
for i = 1:N
    A_t1(4*(i-1)+1:4*i,6*N+6+i) = [1;-1;1;-1];
    A_t1(4*(i-1)+1:4*i,4*(N+1)+1+i) = [g(1)*b1;g(1)*b1;-g(1)*b1;-g(1)*b1];
end
b_t1(1:4*N,1) = b1*wmax(1);
for i = 1:N
    A_t2(4*(i-1)+1:4*i,7*N+6+i) = [1;-1;1;-1];
    A_t2(4*(i-1)+1:4*i,5*(N+1)+1+i) = [g(2)*b2;g(2)*b2;-g(2)*b2;-g(2)*b2];
end
b_t2(1:4*N,1) = b2*wmax(2);

%% Warunek lekkiej zmiany momentu w czasie (t[k+1]-t[k]<g*tmax/N)
A_dt1 = zeros(2*(N-1),10*N+10);
A_dt2 = zeros(2*(N-1),10*N+10);
b_dt1 = zeros(2*(N-1),1);
b_dt2 = zeros(2*(N-1),1);
for i = 1:N-1
    A_dt1(i,6*N+6+i:6*N+6+i+1) = [-1,1];
    A_dt1(N-1+i,6*N+6+i:6*N+6+i+1) = [1,-1];
    A_dt2(i,7*N+6+i:7*N+6+i+1) = [-1,1];
    A_dt2(N-1+i,7*N+6+i:7*N+6+i+1) = [1,-1];
end
b_dt1(1:2*(N-1),1) = g(1)*tmax(1)/(1/2*N);
b_dt2(1:2*(N-1),1) = g(2)*tmax(2)/(1/2*N);
%% Finalna macierz i wektor
A = [A_F;A_t1;A_t2;A_dt1;A_dt2];
b = [b_F;b_t1;b_t2;b_dt1;b_dt2];

end

