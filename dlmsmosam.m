function out=dlmsmosam(dlm, nsample)
%DLMSMOSAM sample from DLM model 
% out=dlmsmosam(dlm, nsample)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<2; nsample = 1; end

[p,m] = size(dlm.F);
n = length(dlm.yhat);
m = m + size(dlm.XX,2)*p;
out = zeros(m,n,nsample);

hw = waitbar(0,'Generating sample');
% there is mcmc chain that has samples from diagonal of W
for isample = 1:nsample
  if fix(isample/10)==isample/10;waitbar(isample/nsample,hw);end
  W = dlm.W;
  V = dlm.V;
  G = dlm.G;
  if isfield(dlm,'chain')
    if isfield(dlm,'vinds') && not(isempty(dlm.vinds))
      vinds = dlm.vinds; % can be only 1 now
    else
      vinds = 0;
    end
    cind = ceil(rand(1)*size(dlm.chain,1));
    if vinds>0
      V = V.* dlm.chain(cind,1);
    end
    for i = 1:length(dlm.winds)
      if dlm.winds(i)>0
        W(i,i) = dlm.chain(cind,dlm.winds(i)+vinds).^2;
      end
    end
    if isfield(dlm,'ginds')
      i0 = max(dlm.winds)+vinds;
      for i=1:length(dlm.ginds)
        G(dlm.ginds(i),dlm.ginds(i))=dlm.chain(cind,i0+i);
      end
    end
  end
  x0 = dlm.x(:,1);
  C0 = dlm.C(:,:,1);
 % x00 = mvnorrnan(1,x0,C0);
  o = dlmsmo(dlm.y,dlm.F,V,x0,G,W,C0,dlm.XX,1,1,0);
  out(:,:,isample) = o.xr;
end

delete(hw);