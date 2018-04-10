function lik = armalik(y,phi,th,s)
%ARMALIK -2*log(likelihood) of ARMA(p,q) model
% lik = armalik(phi,theta,sig)
% phi, AR coefficients
% theta, MA coefficients
% sig, innovation standard deviation

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2014/12/28 $

[G,F,W,C0]=armasys(phi,th);
V = zeros(size(y));
x0 = zeros(size(G,2),1);

if exist('dlmmex') == 3
  % mex boosted likelihood calculations
  lik = dlmmex(y,F,V,x0,G,W*s^2,C0*s^2,[]);
else
  lik = getfield(dlmsmo(y,F,V,x0,G,W*s^2,C0*s^2,[],0,0),'lik');
end
