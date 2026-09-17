load dom_lshape_16.mat
%load dom_lshape_32.mat
G = D1.G;
L = D1.L;
mesh = meshgen_forward_adaptive(D1,0.1,100*D1.L,1);
p = mesh.p; 
t = mesh.t; 
e = mesh.e; %Knoten, Kanten, Dreiecke
figure(1)
pdeplot(p,e,t);
hold on
plot(p(1,:),p(2,:),'x');
counter = 0;
anzahl_elektroden = max(G(:,3));
m = sum(G(:,3) == 1);

for j = 1:length(G)
    if G(j,3) ~= 0     
        counter = counter+1;
        plot(G(j,1),G(j,2),'o'); %plot electroden
        if counter ~= 0 && mod(counter,m) == 0
            plot([G(j,1),G(j-m+1,1)],[G(j,2),G(j-m+1,2)],'color','green','LineWidth',5);
            hold on
        end
    end
end


points = (mesh.p(:,mesh.t(1,:))+mesh.p(:,mesh.t(2,:))+mesh.p(:,mesh.t(3,:)))./3;
xi = [-0.5;-0.5];
cond = zeros(length(points),1);
for j = 1: length(points)
    cond(j,1) = sigma([points(1,j);points(2,j)],xi);
end
%length
E = zeros(anzahl_elektroden,1);
for j = 1:anzahl_elektroden
    elec1 = G(:,3) == j;
    elec11 = G(elec1,1);
    elec12 = G(elec1,2);
    for i = 1:length(elec11)-1
        E(j) = E(j) + sqrt((elec11(i+1)-elec11(i))^2 + (elec12(i+1)-elec12(i))^2);  
    end
end

I1 = -ones(L,1)/L; I1(12) = (L-1)/L;
I2 = zeros(L,1); I2(3) = 1; I2(7) = -1;

I = [I1,I2];

z=0.04*ones(L,1);

I1 = -ones(L,1)/L;
I1(12) = (L-1)/L;
I2 = zeros(L,1); 
I2(3) = 1; 
I2(7) = -1;

[u,U] = solve_CEM(I,cond,E,G,z,mesh);


figure(2);
s1 = subplot(2,2,1);
pdeplot(p,e,t,'XYData',u(:,1), 'ZData',u(:,1), 'Mesh', 'on'); colormap('jet')
hold on

for j = 1:length(G)
    if G(j,3) ~= 0     
        counter = counter+1;
        plot(G(j,1),G(j,2),'o'); %plot electroden
        if counter ~= 0 && mod(counter,m) == 0
            plot([G(j,1),G(j-m+1,1)],[G(j,2),G(j-m+1,2)],'color','green','LineWidth',5);
            hold on
        end
    end
end


caxis([-1.7,1.7]);
s2 = subplot(2,2,2);
pdeplot(p,e,t,'XYData',u(:,2), 'ZData',u(:,2), 'Mesh', 'on'); colormap('jet')
caxis([-1.7,1.7]);
s3 = subplot(2,2,3);
plot(1:anzahl_elektroden,U(:,1));
hold on
plot(1:anzahl_elektroden,I(:,1));
s4 = subplot(2,2,4);
plot(1:anzahl_elektroden,U(:,2));
hold on
plot(1:anzahl_elektroden,I(:,2));
hold off

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


function output = sigma(location, xi)
    if norm(location-xi,2)<0.25
        output = 5;

    elseif norm(location - xi,2)<0.5
        output = 0.05;
    else
        output = 1;  
    end
end

%Spannung misst, Ströme Eingegeben