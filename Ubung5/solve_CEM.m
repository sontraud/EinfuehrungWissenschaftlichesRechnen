function [u,U] = solve_CEM(I,cond,E,D1G,z,mesh)
    %Eingabeparameter:
    %I (LxM)-Matrix für Ströme
    %cond diskrete Leitfähigkeit
    %E Längen der ELektroden
    %D1G Knoten auf Rand des Gebiets + Lage der Elektroden
    %z Kontaktimpedanzen
    %mesh Triangulierung des Gebiets
    %Ausgabe:
    %u (nxM)-Matrix mit n = #Knoten/Schwerpunke für innere Potentiale
    %U (LxM)-Matrix für Elektrodenspannungen

    M = size(I,2); %Anzahl Ströme
    L = size(I,1); %Anzahl Elektroden

    [A1,~,~] = assema(mesh.p,mesh.t,cond,0,0); %Berechne A1 als Steifigkeitsmatrix gemäß (9.21)
    D = sparse(1:L,1:L,E./z,L,L);
    [A2,B] = CEM_FEmatrix(mesh.p,mesh.e,D1G,L,z); %Initialisiere Matrizen A2 und B
    n = size(A1,1); %#Knoten/Schwerpunkte
    Matrix = [A1+A2,B;B',D;sparse(1,n),sparse(ones(1,L))]; %Matrix zum LGS
    
    u = zeros(n,M);
    U = zeros(L,M);

    for i = 1:M
        vec = [sparse(n,1);sparse(I(:,i));sparse(1,1)]; %Rechte Seite des LGS
        vhelp = full(Matrix\vec); %Lösung LGS, anschließend Aufteilen auf u und U
        u(:,i) = vhelp(1:n);
        U(:,i) = vhelp((n+1):(n+L));
    end
end