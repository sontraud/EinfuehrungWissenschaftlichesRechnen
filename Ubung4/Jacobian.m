function Jac = Jacobian(UE,IM,mf)
%Assemblierung der Jacobi-Matrix bez. Leitfähigkeit auf dem FE-Gitter
%
%  UE  innere Potentiale bez. der virtuellen Ströme I_\diamond
%  IM  tatsächlich angelegte Ströme
%  mf  FE-Gitter-Struktur (erzeugt mit meshgen_forward_adaptive.m)
%  
%  Jac ist eine (M*L,length(mf.t)) Matrix, d.h. die Ausgabe für 
%      mehrere Ströme erfolgt untereinander!

[L,M] = size(IM);

UMx = zeros(length(mf.points),M);     %storage for the gradient of the solution
UMy = zeros(length(mf.points),M);     %storage for the gradient of the solution
UEx = zeros(length(mf.points),L);     %storage for the gradient of the solution
UEy = zeros(length(mf.points),L);     %storage for the gradient of the solution

% virtual setting
for j = 1:L
  U1 = UE(mf.t(1,:),j);
  U2 = UE(mf.t(2,:),j);
  U3 = UE(mf.t(3,:),j);
  UEx(:,j) = U1.*mf.g1x+U2.*mf.g2x+U3.*mf.g3x;
  UEy(:,j) = U1.*mf.g1y+U2.*mf.g2y+U3.*mf.g3y;  
end
clear U1 U2 U3 UE;

for k=1:M
  for j=1:L
    UMx(:,k) = UMx(:,k) + IM(j,k)*UEx(:,j);
    UMy(:,k) = UMy(:,k) + IM(j,k)*UEy(:,j);
  end
end

% assemble the Jacobian
H = aux_dotkron(UMx',UEx')+aux_dotkron(UMy',UEy');
clear UMx UMy UEx UEy;

Jac = -aux_dotkron(H,mf.areas');