function [G,F,W,C0] = armasys(phi,th)
%ARMASYS generate state space representation of ARMA model
% [G,F,W,C0] = armasys(phi,theta)

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2014/12/28 $

p = length(phi);
q = length(th);
pq = max(p,q+1);
G = [[phi(:);zeros(pq-p,1)],[eye(pq-1);zeros(1,pq-1)]];
F = [1 zeros(1,pq-1)];
R = [1;th(:);zeros(pq-q-1,1)];
W = R*R';

% C0, initial state covariance
if nargout>3
  C0 = reshape((eye(pq^2)-kron(G,G))\W(:),pq,pq);
end
