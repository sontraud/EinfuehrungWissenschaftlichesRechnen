sigma = 10;
beta = 8/3;
rho = 28;
v0=[1;1;1];
t=[0,50];
ode_f = @(t,x) [sigma*(x(2)-x(1));x(1)*(rho-x(3))-x(2);x(1)*x(2)-beta*x(3)];
rho_rk = 0.99;
h0=0.01;
hmax=0.2;
tol = 2*10^(-7);
[t,y,fail,no_fun] = odesolverRK(ode_f,t,v0,h0,hmax,tol,rho_rk);

h=[];
for j = 1:length(t)-1
    h=[h,t(j+1)-t(j)];
end



figure(1);
figure('units','normalized','outerposition',[0 0 1/2 1])
axis([-30 30 -30 30 -10 50]), grid on
hold on;
f=animatedline('Color', 'black', 'LineWidth',1);
view(3)
%hold off

figure(2);
figure('units','normalized','outerposition',[1/2 0 1/2 1])
axis([0 50 0 hmax])
g=animatedline('Color', 'red', 'LineWidth',1);
%axis([0 t(end) 0 hmax]);
%hold on

for k = 1:length(t)-1
    addpoints(f,y(1,k),y(2,k),y(3,k));
    addpoints(g,t(k),h(k))
    drawnow;
end


global_error = norm([y(1,end),y(2,end)]-[v0(1),v0(2)],inf);
fprintf('number of function evaluation %d', no_fun);
fprintf('\n number of failed iterations %d', fail);
fprintf('\n global error %f', global_error);