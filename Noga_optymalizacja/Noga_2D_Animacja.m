%% Animacja dla nogi 2D (symulacja)
% Stałe:
l1 = 0.25;
l2 = 0.25;

% Animacja: 
clf;
l_1 = plot([0,1],[0,0],'Color',[0 0 0],'LineWidth',3);
l_2 = plot([0,1],[0,0],'Color',[0 0 0],'LineWidth',3);
p_0 = plot(0,1,'bo','Color','red','LineWidth',2);
p_1 = plot(0,1,'bo','Color','red','LineWidth',2);
p_2 = plot(0,1,'bo','Color','red','LineWidth',2);
axis([-1 1 -1 1]);
grid on;
n = size(t);
dt = 0.005;
N_a = round(t(n(1))/dt);
for j = 1:N_a
    [u,i] = min(abs(t-ones(n(1),1)*j*dt));
    delete(l_1);  
    delete(l_2);
    delete(p_0)
    delete(p_1);
    delete(p_2);
    x_1 = l1*cos(alfas(i));
    y_1 = ys(i) + l1*sin(alfas(i));
    x_2 = x_1 + l2*cos(alfas(i)+betas(i));
    y_2 = y_1 + l2*sin(alfas(i)+betas(i));
    hold on;
    l_1 = plot([0,x_1],[ys(i),y_1],'Color',[0 0 0],'LineWidth',3);
    l_2 = plot([x_1,x_2],[y_1,y_2],'Color',[0 0 0],'LineWidth',3);
    p_0 = plot(0,ys(i),'bo','Color','blue','LineWidth',4);
    p_1 = plot(x_1,y_1,'bo','Color','red','LineWidth',2);
    p_2 = plot(x_2,y_2,'bo','Color','red','LineWidth',2);
    pause(dt);
end