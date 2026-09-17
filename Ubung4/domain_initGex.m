function Gex = domain_initGex(G,verbose)
%MESH_INITGEX initialize domain geometry from domain boundary data G

if nargin < 2
  verbose = 0;
end

% this only works if G is parametrized counter-clockwise
seg = size(G,1);

Gex = zeros(7,seg);
Gex(1,:) = 2;
Gex(2,:) = G(:,1)';
Gex(3,:) = Gex(2,[2:end 1]);
Gex(4,:) = G(:,2)';
Gex(5,:) = Gex(4,[2:end 1]);
Gex(6,:) = 1;
Gex(7,:) = 0;

if verbose > 2
  figure;
  pdegplot(Gex,'edgeLabels','on','subdomainLabels','on');
  xlim([min(G(:,1))-0.1 max(G(:,1))+0.1]);
  ylim([min(G(:,2))-0.1 max(G(:,2))+0.1]);
  axis equal;
  title('PDE toolbox boundary geometry');
  pause(0.01);
end

end