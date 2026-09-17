function [p,t,e] = meshgen_refineApprox(Gex,p,t,e,G,nttarget,relLim,verbose)
%MESHGEN_REFINEAPPROX refine approximation mesh

nmax = 20;
n = 0;
ntnew = size(t,2);
bla = tic;

if verbose > 1
  pdemesh(p,e,t);
  pause(1);
end


while n < nmax && ntnew < nttarget
  n = n+1;
  % compute triangle centers
  tcenter = (p(:,t(1,:)) + p(:,t(2,:)) + p(:,t(3,:)))./3;
  
  % compute triangle areas
  [ar,~,~,~,~,~,~] = pdetrg(p,t);
  A0res = 0.1*pi; % arbitrary quantity

  % computed triangle areas according to heuristic
  G3L = logical(G(:,3));
  G3R = logical(G([end 1:end-1],3));
  
  elL = G( G3L & ~G3R,1:2)'; elL = elL(1,:) + 1i*elL(2,:);
  elR = G(~G3L &  G3R,1:2)'; elR = elR(1,:) + 1i*elR(2,:);
  tz = tcenter(1,:) + 1i*tcenter(2,:);

  nt = size(t,2);
  dQP = zeros(1,nt);
  dQM = zeros(1,nt);
  
  
  for k=1:nt
    dQP(k) = min(abs(elR-tz(k)));
    dQM(k) = min(abs(elL-tz(k)));
  end
  
  dQPMmax = max(dQP)*max(dQM);
  rQH = 0.95*1-sqrt((dQP.*dQM)/dQPMmax); % limit radius away from 1
  AQH = 0.5*(A0res + A0res*cos(pi*rQH));
  
  dA = ar./AQH;
  dArel = dA / min(dA);
  
  amin = 1;
  
  ntnew = nt;
  
  ref4 = zeros(nt,1);
  ref2 = zeros(nt,1);
  
  deltan = ntnew+max(6,nttarget-ntnew);
  
  while ntnew < deltan
    [amax,imax] = max(dArel);
    if amax > 6*amin
      ref4(imax) = ref4(imax)+1;
      dArel(imax) = dArel(imax)/4;
      ntnew = ntnew+5;
    else
      break;
    end
  end

  while ntnew < deltan
    [amax,imax] = max(dArel);
    if amax > amin
      ref2(imax) = ref2(imax)+1;
      dArel(imax) = dArel(imax)/2;
      ntnew = ntnew+2;
    else
      break;
    end
  end

  ref4 = (ref4==max(ref4) & ref4>0);
  ref4ind = find(ref4);
  if length(ref4ind)>=1
    % refinemesh is stupid. it interprets a single element as "row vector" (i.e. subdomain to refine)
    if length(ref4ind)==1; ref4ind = [ref4ind; ref4ind]; end
    
    [p,e,t] = refinemesh(Gex,p,e,t,ref4ind,'regular');
    p = jigglemesh(p,e,t,'Opt','minimum','Iter',5);
    ref2ind = [];
  else  
    ref2ind = find(ref2);
    if any(ref2ind)
      [p,e,t] = refinemesh(Gex,p,e,t,ref2ind,'longest');
      p = jigglemesh(p,e,t,'Opt','minimum','Iter',5);
    end
  end
  
  if verbose > 1
    fprintf('ref4=%i, ref2=%i\n',numel(ref4ind), numel(ref2ind));
  end
  
  ntnew = size(t,2);
  
  if verbose > 1
    pdemesh(p,e,t);
    pause(1);
  end
  
end

if verbose
  fprintf('Did %i triangle refinements in %4.3es.\n',n,toc(bla));
end

end






