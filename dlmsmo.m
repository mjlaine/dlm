function out = dlmsmo(y,F,V,x0,G,W,C0,X, sample, smooth, lik)
%DLMSMO DLM smoothing
% out = dlmsmo(y,F,V,x0,G,W,C0, X)
% Uses Kalman smoother for dynamic linear state space model:
%
%  y(t) = F*x(t)   + v
%  x(t) = G*x(t-1) + w,    for t=1:n
%
% with
%
%  x(1) ~ N(x0,C0), v ~ N(0,V), and w ~ N(0,W).
%
% Input:
% y  observations, n*p
% F  observation operator, p*(m-nx)
% V  observation uncertainty std, n*p
% x0 initial state, m*1
% G  system evolution matrix, m*m
% W  model error covariance, m*m
% C0 initial state uncertainty covariance, m*m
% X  external covariate matrix, n*nx (optional)
%
% Output:
% out.x smoothed state m*n
% out.C smoother covariance m*m*n
% etc...

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin < 9
  sample = 1; % generate also one sampled state for MCMC
end

if nargin < 10
  smooth = 1; % calculate smoother, also, not needed for likelihood
end

if nargin < 11
  lik = 1; % calculate the likelihood (and chol(Cp)), can be slow in large models
end

if nargin < 8 || isempty(X)
  X = zeros(size(y,1),0);
end

[p,m] = size(F);  % n_series, n_states
n = size(y,1);    % n_observations
m = m + size(X,2)*p; % covariates, i.e. the proxy variables

x = zeros(m,n); % states
C = zeros(m,m,n); % state uncertainty covariances

% initial values
x(:,1)   = x0;
C(:,:,1) = C0;

% collect these also
Cp = zeros(p,p,n);  % obs prediction uncertainty
v  = zeros(p,n);    % prediction residuals
K  = zeros(m,p,n);  % Kalman gain

% Kalman filter recursion, now x is the one step prediction mean
for i=1:n
  FF = [F,kron(X(i,:),eye(p))];
  ig = not(isnan(y(i,:))); % non NaN y
  v(ig,i) = y(i,ig)' - FF(ig,:)*x(:,i);
  Cp(:,:,i) = FF*C(:,:,i)*FF';
  Cp(ig,ig,i) = Cp(ig,ig,i) + diag(V(i,ig).^2);
%  if any(diag(Cp(:,:,i))<=0), keyboard;end
  Cp(:,:,i)  = triu(Cp(:,:,i)) + triu(Cp(:,:,i),1)'; % fix symmetry ...
  K(:,ig,i) = G*C(:,:,i)*FF(ig,:)'/Cp(ig,ig,i);
  if i<n
    L = G-K(:,:,i)*FF;
    x(:,i+1) = G*x(:,i) + K(:,:,i)*v(:,i);
    C(:,:,i+1) = G*C(:,:,i)*L' + W;
    C(:,:,i+1)  = triu(C(:,:,i+1)) + triu(C(:,:,i+1),1)'; % fix symmetry ...
%    if any(diag(C(:,:,i+1))<=0), disp('C negative'); keyboard;end
  end
end

% save filtered (predicted) state
out.xf = x;
out.Cf = C;

if smooth

  % Smoothing recursion
  r = zeros(m,1);
  N = zeros(m,m);
  if sample
    xr = zeros(m,n);
    yr = zeros(p,n);
  end

  for i=n:-1:1

    FF = [F,kron(X(i,:),eye(p))];

    if sample
      eta = W*r; % mean smoothing disturbance
      Ceta = W-W*N*W;
      xr(:,i) = mvnorrnan(1,eta,Ceta); % disturbance
    end

    L = G-K(:,:,i)*FF;
    ig = not(isnan(y(i,:)));
    FFCp = FF(ig,:)'/Cp(ig,ig,i);
    r = FFCp*v(ig,i) + L'*r;
    N = FFCp*FF(ig,:) + L'*N*L;
    x(:,i) = x(:,i) + C(:,:,i)*r;
    C(:,:,i) = C(:,:,i) - C(:,:,i)*N*C(:,:,i);

    C(:,:,i)  = triu(C(:,:,i)) + triu(C(:,:,i),1)'; % fix symmetry ...
    for ii=1:size(C,1), C(ii,ii,i) = abs(C(ii,ii,i)); end
%    if any(diag(C(:,:,i))<=0), disp('C smo negative'); keyboard;end

  end

  if sample % disturbances to states
    FF = [F,kron(X(1,:),eye(p))];
    out.xrd = xr; % save disturbance sample
    xr(:,1) = mvnorrnan(1,x(:,1),C(:,:,1))';
    yr(:,1) = FF*xr(:,1) + mvnorrnan(1,zeros(1,p),diag(V(1,:).^2))';
    yr(isnan(y(1,:)),1) = NaN;
    for i=2:n
      FF = [F,kron(X(i,:),eye(p))];
      xr(:,i) = G*xr(:,i-1) + xr(:,i);
      ig = not(isnan(y(i,:)));
      yr(~ig,i) = NaN;
      yr(ig,i) = FF(ig,:)*xr(:,i) + mvnorrnan(1,zeros(1,sum(ig)),diag(V(i,ig).^2))';
    end
  end

end % if smooth

yhat = zeros(n,p);
ystd = zeros(n,p);
Csd = zeros(n,m); % smoothing error std
for i=1:n
  FF = [F,kron(X(i,:),eye(p))];
  yhat(i,:) = (FF*out.xf(:,i))'; % filter prediction
  Csd(i,:) = sqrt(diag(squeeze(C(:,:,i))));
  ystd(i,:) = sqrt(diag(FF*C(:,:,i)*FF'+diag(V(i,:).^2)))'; % C is smooth now
end

out.x = x;
out.C = C;
out.xstd = Csd;
out.G = G;
out.F = F;
out.W = W;
out.y = y;
out.V = V;
out.x0 = x0;
out.C0 = C0;
out.XX = X;
out.yhat = yhat; % filter prediction obs
out.ystd = ystd; % smoother prediction obs unc (FIX ME)
out.resid0 = y-yhat; % raw smoother resiaduals
out.resid = out.resid0./V; % scaled smoother residuals
out.ssy = sumnan(out.resid0.^2);
out.v = v'; % filter prediction residual
out.Cp = squeeze(Cp); % filter prediction obs unc
out.s2 = sumnan(out.resid.^2)./(n-sum(isnan(y))); % FIXME

out.nobs = n - sum(isnan(y));
out.lik = NaN;
out.mse = NaN;
out.mape = NaN;

if p==1 % -2*log likelihood for single series case
  out.v(isnan(y)) = NaN;
  if lik
    out.lik = sumnan(out.v.^2./out.Cp + log(out.Cp));
  end
  out.resid2 = v'./sqrt(squeeze(Cp)); % scaled prediction residuals
  out.resid2(isnan(y)) = NaN;
  out.mse = sumnan(out.resid2.^2)./out.nobs;
  out.mape = sumnan(abs(out.resid2)./y)./out.nobs;
else % multivariate calculations
  if lik
    out.lik = 0;
    out.resid2 = NaN(n,p);
    for i=1:n
      ig = not(isnan(y(i,:)));
      Cpchol = chol(Cp(ig,ig,i));
      out.lik = out.lik + v(ig,i)'/Cp(ig,ig,i)*v(ig,i) + sum(log(diag(Cpchol)))*2;
      out.resid2(i,ig) = v(ig,i)'/Cpchol;
    end
    % no resid2 or mse, if lik=0
    out.mse = sumnan(out.resid2.^2)./out.nobs;
    out.mape = sumnan(abs(out.resid2)./y)./out.nobs;
  end
end

if sample && smooth
  % need to run the smoother again to get the sample
  smo = dlmsmo(yr',F,V,x0 ,G,W,C0, X, 0, 1, 0);
  out.xr = xr - smo.x + x;
  out.xrp = xr;
  out.yrp = yr;
  out.ss = sum((out.xr(:,2:end)-G*out.xr(:,1:end-1))'.^2); % state sum-of-squares
end

out.class = 'dlmsmo';
