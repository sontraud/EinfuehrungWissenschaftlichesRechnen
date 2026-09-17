function [A2,B] = CEM_FEmatrix(p,e,G,L,z)
% A2 und B sind Matrizen, die zur Assemblierung der Steifigkeitsmatrix K
% der FEM-Diskretisierung des CEM auftreten, siehe z.B. (3.4)
% in A. Lechleiter und A. Rieder "Newton regularization for impedance tomography:
% a numerical study", Inverse Problems 22, 1967-1987, 2006.
% 
%  K = [A1+A2, B; B', D]
%  A1 und D werden nicht in dieser Routine berechnet.
%
%  Eingabe: p,e Liste der Knoten und (Rand)kanten
%           G  Elektrodenkonfiguration (D1.G aus der Geometriestruktur D1)
%           L Anzahl der Elektroden
%           z Vektor der Kontaktimpedanzen

sC = zeros(size(G(:,3)));        % segment conductivity

for k=1:L
  sC(G(:,3)==k) = 1/z(k);
end


dim=length(p);
idx = 1;              % starting index: first segment for the "ordered" boundary segments in e

% compute boundary terms
segR = zeros(length(e),1);
for j = 1:length(e)
  segR(j) = idx;
  if j < length(e)
    idx = find(e(1,:) == e(2,idx)); % find boundary segments in circular order
  end
end
segL = segR([end 1:end-1]);       % "left-side" segments

pBd  = e(1,segR);                 % points between "left-side" and "right-side" segments

zRHS = sC(e(5,segR)); % inverse conductivity at RHS segment of point
zLHS = sC(e(5,segL)); % inverse conductivity at LHS segment of point

slRHS = sqrt(sum((p(:,pBd([2:end 1])) - p(:,pBd)).^2))';% RHS segment length
slLHS = slRHS([end 1:end-1]);                           % LHS segment length
%clear p;

% ([Polydorides 2002, (3a) 2nd part for i=j+1 and i=j-1])
A2up = (1/6) .* zLHS .* slLHS; %(phi_i and phi_(i-1) only have common support on LHS edge of the point i)
A2dn = (1/6) .* zRHS .* slRHS; %(phi_i and phi_(i+1) only have common support on RHS edge of the point i)

% ([Polydories 2002, (3a) 2nd part for i=j])
A2dg = 2 .* (A2up+A2dn);       %(phi_i^2 has support on LHS and RHS (if they are on an electrode)

% assemble A2
%"upper", "lower" and "diagonal" entry coordinates
A2icoord = [pBd';                pBd';            pBd'];
A2jcoord = [pBd([end 1:end-1])'; pBd([2:end 1])'; pBd'];

A2 = sparse(A2icoord,A2jcoord,[A2up; A2dn; A2dg],dim,dim);
Btmp = 3 .* (A2up+A2dn);


% electrode is stored for each edge segment.
% to find the electrode of a point, check left and right edge segment.
% either they match, or one of them is 0 (not belonging to an electrode)
Points2Elec = max(G(e(5, segL),3),G(e(5, segR),3));

%clear segL segR;

% [Polydorides 2002, (3b)]
Bval = zeros(length(e),L);
[Bicoord,Bjcoord] = meshgrid(1:L,pBd);
for l=1:L
  Bval(:,l) = -Btmp.*(Points2Elec==l);
end
B = sparse(Bicoord(:),Bjcoord(:),Bval(:),L,dim)';
