%i) Generiere Mesh und berechne Schwerpunkte der Dreiecke:
clear all;
load('dom_lshape_16.mat'); %Lese Struktur D1 ein, D1.G enthält Knoten auf dem Rand + Lage der Elektroden
mesh = meshgen_forward_adaptive(D1,0.1,100*D1.L,1); %Generiere Triangulierung
L = size(D1.G,1)/8; %Anzahl Elektroden

mesh.points = (mesh.p(:, mesh.t(1, :)) + mesh.p(:, mesh.t(2, :)) + mesh.p(:, mesh.t(3, :)))./3; %Spaltenweise Schwerpunkte der Dreiecke
location.x = mesh.points(1,:); %location enthält Schwerpunkte
location.y = mesh.points(2,:);
n = length(location.x); %Anzahl Knoten/Schwerpunkte

%ii) Definiere E (Elektrodenlängen), z (Kontaktimpedanzen), c_plus (Exakte Lsg), cond (exakte diskr. Leitfähigkeit
%in Schwerpunkten), I (Ströme), U_exakt (exakte Messwerte), I_virtuell (virtuelle Ströme):
E = zeros(L,1);
for k = 1:L
    for l = 1:3
        E(k) = E(k) + norm(D1.G(8*(k-1)+l+1,1:2) - D1.G(8*(k-1)+l,1:2),2);
    end
end

z = 0.04*ones(L,1);

c_plus = [-0.5;0.5;1];
cond = sigmafunction(location,c_plus);

I = eye(L);
I = I - [I(L,:);I(1:(L-1),:)];
[~,U_exakt] = solve_CEM(I,cond,E,D1.G,z,mesh); %Bestimme exakte Messwerte über exakte Leitfähigkeit
I_virtuell = eye(L) - (1/L)*ones(L,L);

%iii) Bestimme Messwerte mit Rauschen  U_delta:
delta = 0.03; %Rauschpegel
rng(22);
R = randn(L); %Einträge von R sind u.i. (0,1)-normalverteilt
U_delta = U_exakt + delta*(norm(U_exakt,'fro')/norm(R,'fro'))*R;
U_delta_fro = norm(U_delta,'fro'); %Frobeniusnorm von U_delta 

%iv) Definiere c_start und Jacobi-Matrix von Ph_sigma(c):
c_start = zeros(3,1);
Ph_sigma_jac = [sigmafunction(location,[1;0;0]);sigmafunction(location,[0;1;0]);sigmafunction(location,[0;0;1])]' - ones(length(cond),length(c_plus));

%v) Führe Iteration mit Steepest Decent Methode (sd) aus:
k = 0;
ck = c_start; %Startiterierte
e_sd = norm(ck - c_plus,2); %Fehler zum Start 
[u,U_virtuell] = solve_CEM(I_virtuell,sigmafunction(location,ck),E,D1.G,z,mesh); %Innere Potentiale und Elektrodenspannungen bei virtuellen Strömen und Startiterierte
Phi_ck = aux_VEtoVM(U_virtuell,I); %Elektrodenspannungen bei Strömen I und Startiterierte
res_sd = norm(U_delta - Phi_ck,'fro'); %Residuum zum Start
Phi_ck_jac = Jacobian(u,I,mesh)*Ph_sigma_jac; %Berechne Phi'(ck)
sk = (Phi_ck_jac')*reshape(U_delta-Phi_ck,[L^2,1]); %Berechne sk

while res_sd(k+1) > 1.05*delta*U_delta_fro
    k = k+1;
    lambda_k = (norm(sk,2)/norm(Phi_ck_jac*sk,'fro'))^2; %Berechne Faktor lambda_k nach Steepest Decent Methode
    ck = ck + lambda_k*sk; %Berechne neue Iterierte
    [u,U_virtuell] = solve_CEM(I_virtuell,sigmafunction(location,ck),E,D1.G,z,mesh); %Innere Potentiale und Elektrodenspannungen bei virtuellen Strömen und neuer Iterierten
    Phi_ck = aux_VEtoVM(U_virtuell,I); %Elektrodenspannungen bei Strömen I und neuer Iterierten
    e_sd = [e_sd,norm(ck-c_plus,2)]; %Ergänze Fehler von neuer Iterierten
    res_sd = [res_sd,norm(U_delta - Phi_ck,'fro')]; %Ergänze Residuum von neuer Iterierten
    Phi_ck_jac = Jacobian(u,I,mesh)*Ph_sigma_jac; %Berechne neues Phi'(ck)
    sk = (Phi_ck_jac')*reshape(U_delta-Phi_ck,[L^2,1]); %Berechne neues sk
end
c_sd = ck; %Iterierte bei Abbruch

%vi) Führe Iteration mit Minimal Error Methode (minerr) aus:
k = 0;
ck = c_start; %Startiterierte
e_minerr = norm(ck - c_plus,2); %Fehler zum Start
[u,U_virtuell] = solve_CEM(I_virtuell,sigmafunction(location,ck),E,D1.G,z,mesh); %Innere Potentiale und Elektrodenspannungen bei virtuellen Strömen und Startiterierte
Phi_ck = aux_VEtoVM(U_virtuell,I); %Elektrodenspannungen bei Strömen I und Startiterierte
res_minerr = norm(U_delta - Phi_ck,'fro'); %Residuum zum Start
Phi_ck_jac = Jacobian(u,I,mesh)*Ph_sigma_jac; %Berechne Phi'(ck)
sk = (Phi_ck_jac')*reshape(U_delta-Phi_ck,[L^2,1]); %Berechne sk

while res_minerr(k+1) > 1.05*delta*U_delta_fro
    k = k+1;  
    lambda_k = (res_minerr(k)/norm(sk,2))^2; %Berechne Faktor lambda_k nach Minimal Error Methode
    ck = ck + lambda_k*sk; %Berechne neue Iterierte
    [u,U_virtuell] = solve_CEM(I_virtuell,sigmafunction(location,ck),E,D1.G,z,mesh); %Innere Potentiale und Elektrodenspannungen bei virtuellen Strömen und neuer Iterierten
    Phi_ck = aux_VEtoVM(U_virtuell,I); %Elektrodenspannungen bei Strömen I und neuer Iterierten
    e_minerr = [e_minerr,norm(ck-c_plus,2)]; %Ergänze Fehler von neuer Iterierten
    res_minerr = [res_minerr,norm(U_delta - Phi_ck,'fro')];%Ergänze Residuum von neuer Iterierten
    Phi_ck_jac = Jacobian(u,I,mesh)*Ph_sigma_jac; %Berechne neues Phi'(ck)
    sk = (Phi_ck_jac')*reshape(U_delta-Phi_ck,[L^2,1]); %Berechne neues sk
end
c_minerr = ck; %Iterierte bei Abbruch

%vii) Plot der Fehler e_k und der Residuuen r_k für beide Methoden:
figure('Units','normalized','OuterPosition',[0,0.03,1,0.5]);
tiledlayout(1,2);

nexttile;
%subplot(1,2,1)
semilogy(linspace(0,length(e_sd)-1,length(e_sd)),e_sd,'LineWidth',1.5);
hold on;
semilogy(linspace(0,length(e_minerr)-1,length(e_minerr)),e_minerr,'LineWidth',1.5);
title('Fehler e_k bzgl. euklidischer Norm','FontSize',16);
xlabel('Iterationsindex k','FontSize',14);
ylabel('Fehler e_k = ||c^k - c^+||_2','FontSize',14);
legend('Steepest Decent Methode','Minimal Error Methode','location','northeast','FontSize',12);
legend('boxoff');

nexttile;
%subplot(1,2,2)
semilogy(linspace(0,length(res_sd)-1,length(res_sd)),res_sd,'LineWidth',1.5);
hold on;
semilogy(linspace(0,length(res_minerr)-1,length(res_minerr)),res_minerr,'LineWidth',1.5);
title('Residuum r_k bzgl. Frobeniusnorm','FontSize',16);
xlabel('Iterationsindex k','FontSize',14);
ylabel('Residuum r_k = ||U^d^e^l^t^a - Phi(c^k)||_F_r_o','FontSize',14);
legend('Steepest Decent Methode','Minimal Error Methode','location','northeast','FontSize',12);
legend('boxoff');

%viii) Plot der Rekonstruktionsfehler für beide Methoden:
reconerr_sd = sigmafunction(location,c_sd) - cond; %Rekonstruktionsfehler bei Steepest Decent Methode
reconerr_minerr = sigmafunction(location,c_minerr) - cond; %Rekonstruktionsfehler bei Minimal Error Methode

max_reconerr = max(abs([reconerr_sd,reconerr_minerr])); %Betrag des betragsmäßig größten Rekonstruktionsfehlers aus beiden Methoden (nur für Plot)

figure('Units','normalized','OuterPosition',[0,0.45,0.5,0.55]);
str = strcat('Rekonstruktionsfehler zu Steepest Decent Methode mit Fehler e_k_*= ',num2str(e_sd(length(e_sd))),':');
pdeplot(mesh.p,mesh.e,mesh.t,XYData=reconerr_sd,ZData=reconerr_sd,Mesh='on');
axis([-1,1,-1,1,-max_reconerr,max_reconerr]);
colormap jet;
caxis([-max_reconerr max_reconerr]);
title(str,'FontSize',16);
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
axis([-1.25,1.25,-1.25,1.25]);

figure('Units','normalized','OuterPosition',[0.5,0.45,0.5,0.55]);
str = strcat('Rekonstruktionsfehler zu Minimal Error Methode mit Fehler e_k_*= ',num2str(e_minerr(length(e_minerr))),':');
pdeplot(mesh.p,mesh.e,mesh.t,XYData=reconerr_minerr,ZData=reconerr_minerr,Mesh='on');
axis([-1,1,-1,1,-max_reconerr,max_reconerr]);
colormap jet;
caxis([-max_reconerr max_reconerr]);
title(str,'FontSize',16);
xlabel('x_1','FontSize',14);
ylabel('x_2','FontSize',14);
axis([-1.25,1.25,-1.25,1.25]);

%ix) Funktion für Leitfähigkeit sigma in Abhängigkeit von c:
function sigma = sigmafunction(location,c)
    %Übergabe mehrerer Punkte in location möglich, Rückgabe der Funktionswerte als Zeilenvektor
    num = length(location.x);
    middlepoints = [-0.5,-0.5,0.5;0.5,-0.5,-0.5]; %Mittelpunkte für Anzatz von sigma
    radians = [0.2,0.25,0.3]; %Radien für Ansatz von sigma
    sigma = ones(1,num);
    for i = 1:3
        sigma(sqrt((location.x(:) - middlepoints(1,i)).^2 + (location.y(:) - middlepoints(2,i)).^2)<radians(i)) = 1 + c(i);
    end
end