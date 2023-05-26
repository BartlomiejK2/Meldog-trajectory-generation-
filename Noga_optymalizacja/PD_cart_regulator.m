%% Funkcja generująca moment napędowy z regulatora kartezjańskiego PD
function [t] = PD_cart_regulator(q_r,q,dq_r,dq,l) % q_r, dq_r - współrzędne i prędkości złączowe referencyjne, q, dq - współrzędne i prędkości złączowe rzeczywiste
k_1 = 200;  %200;
k_2 = 200; %200;
d_1 = 20; %20;
d_2 = 120; %150;
K = [k_1,0;
    0,k_2];
D = [d_1,0;
    0,d_2];
r_r = forward_kin_wahadlo(q_r,l); % położenie końcówki referencyjne względem członu 1 
r = forward_kin_wahadlo(q,l); % położenie końcówki rzeczywiste względem członu 1 
dr_r = forward_vel_wahadlo(q_r,dq_r,l); % prędkość końcówki referencyjna względem członu 1 
dr = forward_vel_wahadlo(q,dq,l); % prędkość końcówki rzeczywista względem członu 1 
J = Jacobi_wahadlo(q,l); % Jakobian względem członu 1
t = J'*(K*(r_r-r)+D*(dr_r-dr));
end

