function y=mvnorrnan(n,mu,C,R)
% Multivariate normal random numbers with possible zero variances.
% C can be singular, also

% marko.laine@fmi.fi, 2013

if size(C,1)~=length(mu)
  error('sizes of mu and C do not match')
end

% first remove zero variances
inds = diag(C) > 2*eps;

y = repmat(mu(:)',n,1);

C1 = C(inds,inds);
mu1 = mu(inds);
n1 = length(mu1);

%[R,q]  = chol(C1); % first try chol
q=1;
if q>0 % then svd
 % warning('sinqular cov');
  [u,s,~] = svd(C1);
  k = find(diag(s)>1e-10,1,'last');
  R = (u(:,1:k)*sqrt(s(1:k,1:k)))';
  n1 = k;
end

%y(:,inds) = bsxfun(@plus,randn(n,n1)*R,mu1(:)');
y(:,inds) = y(:,inds) + randn(n,n1)*R;
