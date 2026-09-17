function VM = aux_VEtoVM(VE,IM)
%AUX_VETOVM get actual potentials from virtual potentials by linear
%combination

[L,M] = size(IM);

VM = zeros(size(VE,1),M);
for k=1:M
  for j=1:L
    VM(:,k) = VM(:,k)+IM(j,k)*VE(:,j);
  end
end

% should be normalized, but might have cumulated round-off errors
VM = VM - (1/L)*repmat(sum(VM),L,1);
end

