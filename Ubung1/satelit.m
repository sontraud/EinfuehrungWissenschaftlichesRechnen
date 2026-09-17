clf
mu = 0.012277471; %Verhältnis Mondmasse zu Gesamtsystem
mu_hut = 1- mu;
t=[0,17.0652166]; %2*t(2) für doppelte Umlaufzeit
h0=2.5*10^(-4); %Startschrittweite
hmax=0.2; %Maximale Schrittweite
tol = 10^(-5); %Genauigkeitsparameter
rho=0.99; %Sicherheitsfaktor in (0,1)
v0=[0.994,0,0,-2.001585106]'; %Startposition, Startgeschwindigkeit
dist_mond = 1/(1+mu); %dist_erde+dist-mond = 1, auf x_1 Achse
dist_erde = 1-1/(1+mu); %dist_erde*masse_erde = dist_mond*masse_mond, %mu = masse_mond/masse_erde +masse_mond
ode_f = @(t,x) [x(3);x(4); x(1)+2*x(4)-mu_hut*(x(1)+mu)/(((x(1)+mu)^2+x(2)^2)^(3/2))-mu*(x(1)-mu_hut)/((x(1)-mu_hut)^2+x(2)^2)^(3/2); x(2)-2*x(3)-mu_hut*(x(2)/((x(1)+mu)^2+x(2)^2)^(3/2))-mu*x(2)/((x(1)-mu_hut)^2+x(2)^2)^(3/2)];
%x_1 Ort Koordinate 1, x_2 Ort Kordinate 2, x
% _3 Ableitung in x_1, x_4
%Ableitung in x_2


[t,y,fail,no_fun] = odesolverRK(ode_f,t,v0,h0,hmax,tol,rho);
h=[];
for j = 1:length(t)-1
    h=[h,t(j+1)-t(j)]; %Schrittweite gegen Zeitschrittlänge
end

h1 = figure(1);
figure('units','normalized','outerposition',[0 0 1/2 1])
plot([dist_mond,dist_erde],[0,0],'o','MarkerFaceColor','green','MarkerEdgeColor','blue');
axis([-1.5 1 -1.5 1.5]), grid on
hold on;
f=animatedline('Color', 'black', 'LineWidth',1);
%hold off

h2 = figure(2);
figure('units','normalized','outerposition',[1/2 0 1/2 1])
g=animatedline('Color', 'red', 'LineWidth',1);
axis([0 t(end) 0 hmax]);
%hold on

for k = 1:length(t)-1
    addpoints(f,y(1,k),y(2,k));
    addpoints(g,t(k),h(k))
    drawnow;
end

global_error = norm([y(1,end),y(2,end)]-[v0(1),v0(2)],inf);
fprintf('number of function evaluation %d', no_fun);
fprintf('\n number of failed iterations %d', fail);
fprintf('\n global error %f', global_error);

%animatedline(y(1,:),y(2,:));