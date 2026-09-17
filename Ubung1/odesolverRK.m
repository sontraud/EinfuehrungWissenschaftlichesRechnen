function [t,y,fail,no_fun] = odesolverRK(ode_f, tspan, v0,h0,hmax,tol,rho)
    k_1 = @(x,y,h) ode_f(x,y); %Tableau für RK4(3)
    k_2 = @(x,y,h) ode_f(x+0.5*h, y+h*1/2*k_1(x,y,h));
    k_3 = @(x,y,h) ode_f(x+0.5*h,y+h*0.5*k_2(x,y,h));
    k_4 = @(x,y,h) ode_f(x+h,y+h*k_3(x,y,h));
    k_5 = @(x,y,h) ode_f(x+h,y+h*(1/6*k_1(x,y,h)+1/3*k_2(x,y,h)+1/3*k_3(x,y,h)+1/6*k_4(x,y,h)));

    %Initialisierung
    x=tspan(1); %Startzeit
    y=v0;
    h=h0;
    b=tspan(2); %Endzeit
    i=1;
    t_ausgewertet = [tspan(1)];
    y_ausgewertet = [v0];
    fail = 0; %verworfene Schrittweitenvorschläge
    no_fun = 1; %Funktionsauswertungen
    k_1_ausgewertet = k_1(x,y,h); %Die erste FUnktionsauswertung
    while x<b
        k_2_ausgewertet = k_2(x,y,h);
        k_3_ausgewertet = k_3(x,y,h);
        k_4_ausgewertet = k_4(x,y,h);
        k_5_ausgewertet = k_5(x,y,h); %k_1 muss nicht ausgewertet werden da Fehlberg Trick, (= altes k_5)

        no_fun = no_fun+4;

        %y=y+h*(1/6*k_1_ausgewertet+1/3*k_2_ausgewertet+1/3*k_3_ausgewertet+1/6*k_4_ausgewertet);
        y_hut = y+h*(1/6*k_1_ausgewertet+1/3*k_2_ausgewertet+1/3*k_3_ausgewertet+1/6*k_4_ausgewertet); %Verfahren mit höherer Ordnung

        fehler = norm(h*1/6*k_4_ausgewertet-h*1/6*k_5_ausgewertet,'inf');
        h_quer = min([2*h;hmax;rho*h*(tol/fehler)^(1/4)]);
        if fehler <= tol
            t_ausgewertet=[t_ausgewertet,x+h];
            y_ausgewertet = [y_ausgewertet, y_hut];
            x = x+h;
            y = y_hut;
            h = min(h_quer,b-x);
            i=i+1;
            k_1_ausgewertet=k_5_ausgewertet;
        else
            h = h_quer;
            fail = fail+1;
        end
    end
    t=t_ausgewertet;
    y=y_ausgewertet;
end


% [t,y,fail,no_fun]=odesolverRK(ode_f,tspan,v0,h0,hmax,tol,rho)
% Die Bedeutung der Eingabeparameter ist
% ode_f function handle f  ̈ur die rechte Seite der DGL,
% tspan Zeitintervall als Vektor [tspan(1) tspan(2)],
% v0 Anfangswert zum Zeitpunkt tspan(1),
% h0 Anfangsschrittweite,
% hmax maximale Schrittweite,
% tol Genauigkeitsparameter f  ̈ur die Schrittweitenkontrolle,
% rho Sicherheitsfaktor f  ̈ur die Schrittweitenkontrolle.
% Ausgegeben werden sollen
% t Vektor der verwendeten Zeitpunkte,
% y L  ̈osungsmatrix (Spalten entsprechen den Zeitpunkten in t),
% fail Anzahl der verworfenen Schrittweitenvorschl  ̈age,
% no_fun Anzahl der Funktionsauswertungen.