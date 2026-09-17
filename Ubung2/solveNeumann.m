function [u_h,err]= solveNeumann(n,w,plot_result)
    h=1/(n+1);
    %Funktionen für Poisson Problem und Neumann Randbedingungen
    %-laplace u = f innen, Ableitung nach normalen u = g
    u=@(x,y,w) sin(2*pi*w*((x.^2)+y))+y;
    f=@(x,y,w) sin(2*pi*w*((x.^2)+y)).*(16*(x.^2)*w^2*pi^2+4*pi^2*w^2)-4*pi*w*cos(2*pi*w*((x.^2)+y));
    u_x = @(x,y,w) cos(2*pi*w*((x.^2)+y))*4*x.*pi*w;
    u_y = @(x,y,w) cos(2*pi*w*((x.^2)+y))*2*pi*w+1;
    %erstelle T und T Schlange
    T=sparse(1:n+2,1:n+2,4*ones(1,n+2))+sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
    T(1,2)=-2;
    T(n+2,n+1) = -2;
    T_schlange = sparse(2:n+2,1:n+1,-ones(1,n+1),n+2,n+2)+sparse(1:n+1,2:n+2,-ones(1,n+1),n+2,n+2);
    T_schlange(1,2) = -2;
    T_schlange(n+2,n+1) = -2;
    kron_1 = kron(speye(n+2),T);
    kron_2 = kron(T_schlange,speye(n+2));

    %Erstelle A zur Diskretisierung
    A = 1/h^2*(kron_1+kron_2);
    x=linspace(0,1,n+2);
    y = x';
    [X,Y] = meshgrid(x,y);
    F=f(X,Y,w);
    
    F(1,:) = F(1,:)-2*u_y(x,0,w)/h;
    F(n+2,:) = F(n+2,:)+2*u_y(x,1,w)/h;
    F(:,1) = F(:,1)-2*u_x(0,y,w)/h;
    F(:,n+2) = F(:,n+2)+2*u_x(1,y,w)/h;

    %zur Projektion ins Bild von A_h
    d=ones(n+2,n+2);
    d(1,:) = 1/2;
    d(n+2,:)=1/2;
    d(:,1) = 1/2;
    d(:,n+2) = 1/2;
    d(1,1)=1/4;
    d(1,n+2)=1/4;
    d(n+2,1)=1/4;
    d(n+2,n+2) = 1/4;
    %Mittlelwert 0
    F=F-((sum(sum(F.*d)))/sum(sum(d.*d)))*d;


    rhs=reshape(F',[(n+2)^2 1]);
    F=[rhs;0];

    A=vertcat(A,ones(1,(n+2)^2));
    u_h = A\F; 
    u_exakt = u(X,Y,w);
    u_exakt = u_exakt - sum(sum(u_exakt))/((n+2)^2)*ones(n+2,n+2);
    err=norm(u_h-reshape(u_exakt',[(n+2)^2 1]),'inf')
    if strcmp(plot_result, 'true')
        u_h_plot = reshape(u_h,[n+2 n+2])';
        %Plotte Lösung
        f1=figure(1);
        f1.Units = 'normalized';
        f1.OuterPosition = [0 0 1/2 1];
        surf(X,Y,u_h_plot);
        title('Lösung u');
        
        %figure(1,'units','normalized','outerposition',[0 0 1/2 1])
        %Plotte Fehler
        surf(X,Y,u_h_plot);
        f2=figure(2);
        f2.Units = 'normalized';
        f2.OuterPosition = [1/2 0 1/2 1];

        surf(X,Y,-u_exakt+u_h_plot);
        title('Fehler')
    end
end