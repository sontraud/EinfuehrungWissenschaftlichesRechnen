function mf = meshgen_forward_adaptive(DOM,hmax,nttarget,verbose,Hgrad)
%MESHGEN_FORWARD_ADAPTIVE creates and adaptively refined triangulation of
%the domain

G = DOM.G;
L = DOM.L;
relLim = 1-2*min(2*DOM.omega)./DOM.circumference; % smallest relative electrode size

Gex = domain_initGex(G,0);
if nargin < 5
  Hgrad = 1.2;
end

if hmax > 0
  [p,e,t] = initmesh(Gex,'MesherVersion','R2013a','Hgrad',Hgrad,'Hmax',hmax);
else
  [p,e,t] = initmesh(Gex,'MesherVersion','R2013a','Hgrad',Hgrad);
end

[p,t,e] = meshgen_refineApprox(Gex,p,t,e,G,nttarget,relLim,verbose);

[ar,g1x,g1y,g2x,g2y,g3x,g3y]=pdetrg(p,t);

mf.p = p; mf.e = e; mf.t = t; clear p e t;
mf.areas = ar'; 
mf.g1x = g1x'; mf.g2x = g2x'; mf.g3x = g3x'; clear g1x g2x g3x ar;
mf.g1y = g1y'; mf.g2y = g2y'; mf.g3y = g3y'; clear g1y g2y g3y;

if verbose
  fprintf('Total number of triangles: %i (target = %i)\n',size(mf.t',1),nttarget);

end


