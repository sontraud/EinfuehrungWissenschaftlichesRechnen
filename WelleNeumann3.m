close all;
t_start = 0;
tmax=40;
tau = 1/20;
n=99;
%n=29;
h=1/(n+1);


m = 800; %tau = (t_end-t_start)/m
tstart = 0;
%t=linspace(t_start,tmax,m);


% u_alt=@(x,y,w) sin(2*pi*w*((x.^2)+y))+y;
% u=@(t,x,y,w) cos(4*pi*t)*u_alt(x,y,w);

f=@(t,x,y) 300*exp(-10*t)*exp(-200*((x-1/2).^2+(y-1/2).^2));
% u_x = @(t,x,y,w) cos(4*pi*t)*cos(2*pi*w*((x.^2)+y))*4*x.*pi*w;
% u_y = @(t,x,y,w) cos(4*pi*t)*(cos(2*pi*w*((x.^2)+y))*2*pi*w+1);
% u_0 = @(t,x,y,w) u_alt(x,y,w);
%u_t_1 =@(t,x,y,w) 0;

T=sparse(1:n+2,1:n+2,4*ones(1,n+2))+sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
T(1,2)=-2;
T(n+2,n+1) = -2;
T_schlange = sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
T_schlange(1,2) = -2;
T_schlange(n+2,n+1) = -2;
kron_1 = kron(speye(n+2),T);
kron_2 = kron(T_schlange,speye(n+2));
A = 1/h^2*(kron_1+kron_2);


%h=animatedline;
x=linspace(0,1,n+2);
y = x';
[X,Y] = meshgrid(x,y);
F = @(t) f(t,X,Y);


%x=A\B löst Ax = B daher B*A^-1

M_plus = sparse(1:(n+2)^2,1:(n+2)^2,ones(1,(n+2)^2)) + tau^2/4*A;
M_minus = sparse(1:(n+2)^2,1:(n+2)^2,ones(1,(n+2)^2))-tau^2/4*A;
[L,U,P,Q] = lu(M_plus); 


%q_alt = reshape(u(0,X,Y,w)',[(n+2)^2 1]);
q_alt = zeros((n+2)^2,1);
p_alt = zeros((n+2)^2,1);
set(gca,'nextplot','replacechildren');
myVideo = VideoWriter("WelleNeumann3.avi");
open(myVideo);
% h2 = figure(2);
% figure('units','normalized','outerposition',[1/2 0 1/2 1])
% g=animatedline('Color', 'red', 'LineWidth',1);
% axis([0 tmax 0 0.05]);
F_neu = F(0);

F_neu(1,:) = 0;
F_neu(n+2,:) = 0;
F_neu(:,1) = 0;
F_neu(:,n+2) = 0;
for j = 0:m-1
    t=j*tau;
    F_akt = F_neu;
    % F_akt = F(t,w);
    F_neu = F(t+tau);
    % 
    % F_akt(1,:) = F_akt(1,:)-2*u_y(t,x,0,w)/h;
    % F_akt(n+2,:) = F_akt(n+2,:)+2*u_y(t,x,1,w)/h;
    % F_akt(:,1) = F_akt(:,1)-2*u_x(t,0,y,w)/h;
    % F_akt(:,n+2) = F_akt(:,n+2)+2*u_x(t,1,y,w)/h;
    
    F_neu(1,:) = 0;
    F_neu(n+2,:) = 0;
    F_neu(:,1) = 0;
    F_neu(:,n+2) = 0;
    
    F_schlange = reshape(F_akt' + F_neu',[(n+2)^2,1]);
    %zw = L\(P*(tau*A*q_alt + M*p_alt + tau/2*(F_schlange)));
    %p_neu = Q*(U\zw);

    %[L,U,P,Q] = lu(S); P*S*Q=L*U; S*q_alt_mod = q_neu
    %((n+2)^2)+tau^2/4*A)*q_neu = q_alt_mod
    %y_neu = Q*U^-1*L^-1*P*q_alt_mod
    zw = L\(P*(M_minus*q_alt + tau*p_alt+tau^2/4*F_schlange));
    q_neu = Q*(U\zw);
    zw = L\(P*(-tau*A*q_alt+M_minus*p_alt+tau/2*F_schlange));
    p_neu = Q*(U\zw);

    %zw = L\(P*(M*q_alt + tau*p_alt+tau^2/4*F_schlange));
    %q_neu = Q*(U\zw);
    
    
    f1=figure(1);
    f1.Units = 'normalized';
    f1.OuterPosition = [0 0 1/2 1];

    surf(X,Y,reshape(q_alt',[n+2,n+2]),'EdgeColor','none');
    %axis([0 1 0 1  5]);
    title(j)
    drawnow
    frame=getframe(gcf);
    writeVideo(myVideo,frame);
    %surf(X,Y,reshape(q_alt,[n+2, n+2]));
    % f2=figure(2);
    % f2.Units = 'normalized';
    % f2.OuterPosition = [1/2 0 1/2 1];
    %addpoints(g,t,norm(reshape(f(t,X,Y,w)',[(n+2)^2,1])-q_alt,'inf'))
    %drawnow;

    p_alt = p_neu;
    q_alt = q_neu;

    %u_exakt = u(t,X,Y,w);
    % plot(t,norm(reshape(u_exakt,[(n+2)^2,1])-q_alt,'inf'));
end
surf(X,Y,reshape(q_alt',[n+2,n+2]),'EdgeColor','none');
    %axis([0 1 0 1 -5 5]);
    title(j)
    drawnow
    frame=getframe(gcf);
    writeVideo(myVideo,frame);
% addpoints(g,t,norm(reshape(u(t,X,Y,w)',[(n+2)^2,1])-q_alt,'inf'))
% drawnow;
close(myVideo);


%surf(X,Y,reshape(q_alt,[n+2, n+2]));
%nicht eindeutige Lösung da nur über Ableitungen