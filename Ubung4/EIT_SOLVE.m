load dom_lshape_16.mat;
G = D1.G;
L = D1.L;
mesh = meshgen_forward_adaptive(D1,0.1,100*D1.L,1);
p = mesh.p; 
t = mesh.t; 
e = mesh.e; %Knoten, Kanten, Dreiecke

c_plus = [-0.5;0.5;1];

points = (mesh.p(:, mesh.t(1, :)) + mesh.p(:, mesh.t(2, :)) + mesh.p(:, mesh.t(3, :)))./3;
mesh.points= points;
anzahl_elektroden = max(G(:,3));

I = zeros(anzahl_elektroden,anzahl_elektroden);
I(1,anzahl_elektroden) = -1;
I(anzahl_elektroden,anzahl_elektroden) = 1; 
for i = 1:anzahl_elektroden-1
    I(i,i) = 1;
    I(i+1,i) = -1;
end
IV = I - (1\length(I));


c = [-0.5; 0.5;1];

z=0.04*ones(L,1);

E = zeros(anzahl_elektroden,1);
for j = 1:anzahl_elektroden
    elec1 = G(:,3) == j;
    elec11 = G(elec1,1);
    elec12 = G(elec1,2);
    for i = 1:length(elec11)-1
        E(j) = E(j) + sqrt((elec11(i+1)-elec11(i))^2 + (elec12(i+1)-elec12(i))^2);  
    end
end
% 
% cond = zeros(length(points),1);
% for j = 1: length(points)
%     cond(j,1) = sigma([points(1,j);points(2,j)],c(1),c(2),c(3));
% end

[u_V,~] = solve_CEM(IV,cond_vec',E,G,z,mesh); %exakte Lösung
cond_vec = sigma_vec(points,c_plus(1),c_plus(2),c_plus(3));
% figure(2)
% pdeplot(mesh.p,mesh.e,mesh.t,XYData=cond_vec,ZData=cond_vec,Mesh= 'on')
%plot(points(1,:),points(2,:),cond_vec);

[u,U] = solve_CEM(I,cond_vec',E,G,z,mesh); %exakte Lösung
delta = 0.3;
rng(22);
R = randn(anzahl_elektroden,anzahl_elektroden);
U_delta = U+ delta*norm(U,'fro')*R/norm(R,'fro');
Jac_ph = [sigma_vec(points,1,0,0)',sigma_vec(points,0,1,0)',sigma_vec(points,0,0,1)'];
zw = size(Jac_ph);
Jac_ph = Jac_ph - ones(zw(1),zw(2));
P_h = sigma_vec(points,c(1),c(2),c(3));
zw = size(P_h);
P_h = 3*ones(zw(1),zw(2))+P_h;

F_I_strich = Jacobian(u_V,I,mesh);
delta = 0.03;



% for i = 1:L
%     plot(1:L,U(:,i)); hold on;
% end
c_k = [0;0;0];
e1 =  [];
counter = 0;
%Phi_k = U*(sigma_vec(points,c_k(1),c_k(2),c_k(3)));
while norm(U_delta - U,'fro') > 1.05*delta*norm(U_delta, 'fro')
    %disp('called')
    cond_vec = sigma_vec(points,c_k(1),c_k(2),c_k(3));
    [u_V,~] = solve_CEM(IV,cond_vec',E,G,z,mesh); %exakte Lösung
    [u,U] = solve_CEM(I,cond_vec',E,G,z,mesh);
    %U_delta = U+ delta*norm(U,'fro')*R/norm(R,'fro');
    % P_h = sigma_vec(points,c_k(1),c_k(2),c_k(3));
    % zw = size(P_h);
    % P_h = 3*ones(zw(1),zw(2))+P_h;
    F_I_strich = Jacobian(u_V,I,mesh);
    Phi_strich = F_I_strich*Jac_ph;
    
    %Phi_k = U*sigma(points,c_k(1),c_k(2),c_k(3));
    %s_k = (Jacobian(u,I,mesh)*Jac_ph)'*reshape((U_delta-Phi_k),256,1);
    s_k = (F_I_strich*Jac_ph)'*reshape((U_delta-U),256,1);
    lambda_k = (norm(s_k,2)^2)./(norm(Phi_strich*s_k,'fro').^2);
    c_k = c_k + lambda_k*s_k;
    e_k = norm(c_k - c_plus,2);
    e1 = [e1, e_k];
    counter = counter+1;
    disp(norm(U_delta-U,'fro'));
end
figure(3)
semilogy(1:length(e1),e1)




function [u,U] = solve_CEM(I,cond,E,G,z,mesh)
    [L,M] = size(I);
    [A1,~,~] = assema(mesh.p,mesh.t,cond',0,1);

    [A2,B] = CEM_FEmatrix(mesh.p,mesh.e,G,L,z);
    D = E./z;
    D = sparse(diag(D));
    A = A1+A2;
    n = size(A1,1);
    S = [A, B ; transpose(B),D; sparse(1,n),sparse(ones(1,L))];
    
    for i =1:M
        vec = [sparse(n,1);sparse(I(:,i));sparse(1,1)];
        vhilfe = full(S\vec);
        u(:,i) = vhilfe(1:n);
        U(:,i) = vhilfe((n+1):(n+L));
    end
end

function c = sigma_vec(x,c1,c2,c3)
    c = ones(1,length(x));
    c(sqrt((x(1,:) +0.5).^2 +(x(2,:) -0.5).^2)<=0.2) = 1+c1; %Kreise disjunkt
    c(sqrt((x(1,:) +0.5).^2 +(x(2,:) +0.5).^2)<=0.2) = 1+c2;
    c(sqrt((x(1,:) -0.5).^2 +(x(2,:) +0.5).^2)<=0.2) = 1+c3;
end
% function output = sigma(x,c1,c2,c3)
%     output = 1;
%     c = ones(1,length(x));
%     %c(sqrt(x(:,1) +0.5)^2 +sqrt(x(:,2) +0.5)<=0.2) = 1+c1;
%     if abs(x-[-0.5;0.5])<=0.2
%         output = output + c1;
%     elseif abs(x-[-0.5;-0.5]) <=0.25
%         output = output + c2;
%     elseif abs(x-[0.5;-0.5]) <=0.3
%         output = output+c3;
%     end
% end