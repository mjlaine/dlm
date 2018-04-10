function out=dlmtsfit(t,y,s,X,options)
% Fit a time series model

% t time in matlab format
% y n*p data matrix
% s n*p std uncertainty
% X n*nx proxy variables

if nargin < 4
  X = [];
end
if nargin < 5
  options = [];
end

% remove NaN's from the begining and end
i0 = isnan(s) | s<=0;
y(i0) = NaN;
if size(y,2)>1
  i1 = find(not(all(isnan(y'))'),1,'first')-1;
  i2 = find(not(all(isnan(y'))'),1,'last')+1;
else
  i1 = find(not(isnan(y)),1,'first')-1;
  i2 = find(not(isnan(y)),1,'last')+1;
end
y([1:i1,i2:end],:) = [];
s([1:i1,i2:end],:) = [];
t([1:i1,i2:end],:) = [];

% scale y with a global scale for all columns
ys = meannan(stdnan(y));
y = y./ys;
s = s./ys;

o = options;
o.trig  = getopt(o,'trig',1);
o.ns    = getopt(o,'ns',12);
o.order = getopt(o,'order',1);
o.arphi = getopt(o,'arphi',[]);
o.trend = getopt(o,'trend',struct());
o.seas = getopt(o,'seas',struct());
o.ar = getopt(o,'ar',struct());

o.level = getopt(o,'level',struct());

o.level.sig  = getopt(o.level,'sig',0);
o.level.clen = getopt(o.level,'clen',0);
o.level.tau = getopt(o.level,'tau',0);

o.trend.sig  = getopt(o.trend,'sig',0.001);
o.trend.clen = getopt(o.trend,'clen',5);
o.trend.tau = getopt(o.trend,'tau',0);
o.seas.sig   = getopt(o.seas,'sig',0.001);
o.seas.clen  = getopt(o.seas,'clen',5);
o.seas.tau  = getopt(o.seas,'tau',0);
o.ar.sig     = getopt(o.ar,'sig',0.001);
o.ar.clen    = getopt(o.ar,'clen',5);
o.ar.tau    = getopt(o.ar,'tau',0);

o.solar = getopt(o,'solar',0);
o.qbo = getopt(o,'qbo',0);

o.refit = getopt(o,'refit',0);


if o.solar
  X = solarfun(t,1);
  % fix missing solar, at the end, hopefully
  ii = find(isnan(X));
  y(ii,:) = [];
  s(ii,:) = [];
  t(ii,:) = [];
  X(ii,:) = [];
end

if o.qbo
  X = [X,qbofun(t,1)];
  ii = find(isnan(X(:,2)));
  y(ii,:) = [];
  s(ii,:) = [];
  t(ii,:) = [];
  X(ii,:) = [];

end

% not used
%o.opt = 0;
%o.fitv = 0; % 
%o.maxfuneval = 1000;
%o.winds =  [0, 1, 2*ones(1,o.trig*2), ones(1,length(o.arphi)>0)*3];
%o.fitar = length(o.arphi)>0; %

% generate system matrices
[G0,F0] = dlmgensys(o);
p = size(y,2);
G = kron(G0,eye(p));
F = kron(F0,eye(p));
n = size(y,1);
q = size(F,2);
q0 = size(F0,2);
nx = size(X,2);


% distance based covariance matrix needed for multicolumn data,
% correlation beween columns
distmat = toeplitz(0:p-1);
cfun = @(d,phi) exp(-abs(d)/phi);
%cfun = @(d,phi) exp(-0.5.*(d/phi).^2);
cmatfun = @(sig,phi,tau) sig.^2*cfun2(distmat,phi) + eye(size(distmat))*tau.^2;

% generate full matrices for all columns
W1 = zeros(p,p);             % level
W1 = cmatfun(o.level.sig,o.level.clen,o.level.tau); % trend
W2 = cmatfun(o.trend.sig,o.trend.clen,o.trend.tau); % trend
W3 = cmatfun(o.seas.sig,o.seas.clen, o.seas.tau); % seas
stack = @(a,b)[a,zeros(size(a,1),size(b,2));zeros(size(b,1),size(a,2)),b];
W = stack(W1,W2); 
for i=1:o.trig
  W = stack(W,W3);
  W = stack(W,W3);
end
if length(o.arphi)>0          % AR
  W4 = cmatfun(o.ar.sig,o.ar.clen,o.ar.tau);
  W = stack(W,W4);
end

% for the proxies
for i=1:nx
  W = stack(W,W1); % add W
  G(q+1:q+p,q+1:q+p) = eye(p); % grow G, F is taken care of by dlmsmo (fix this)
  q = q + p; % p more states needed
end

% initial values
x0 = zeros(q,1);
x0(1:p) = y(1,:);
x0(find(isnan(y(1,:)))) = 0;
C0 = eye(q,q);

% run dlmsmo
out = dlmsmo(y,F,s,x0,G,W,C0, X, 0, 1, 0);
% again, with new initial values, as in dlmfit
if o.refit
  x0 = out.x(:,1);
  C0 = 1*squeeze(out.C(:,:,1));
  out = dlmsmo(y,F,s,x0,G,W,C0, X, 0, 1, 0);
end

% save some extra elements
if isfield(o,'label')
  out.label = o.label;
end

out.ys = ys;
out.y = y;
out.s = s;
out.time = t;
out.options = o;

p = size(out.F,1);

% find out the indexes of various state elements, level, seas, ar, X
if p==1
  ii = 1;
  out.inds.level = ii;
  ii = ii + o.order + 1;
  out.inds.seas = ii:ii+2*o.trig-1;
  ii = ii+2*o.trig;
  out.inds.ar = ones(length(o.arphi),1)+ii-1;
  ii = length(o.arphi) + ii;
  if size(X,2)>0
    out.inds.X = ii:ii+size(X,2)-1;
  else
    out.inds.X = [];
  end
else

  ii = 1;
  out.inds.level = reshape(ii:p,p,1);
  ii = ii*p + o.order*p + 1;
  out.inds.seas = ii:ii+(2*o.trig)*p-1;
  out.inds.seas = reshape(out.inds.seas,p,numel(out.inds.seas)/p);
  ii = ii+(2*o.trig)*p;
%  out.inds.ar = ones(1,length(o.arphi)*p)+ii-1;
  out.inds.ar = ii:(length(o.arphi)*p)+ii-1;
  out.inds.ar = reshape(out.inds.ar,p,numel(out.inds.ar)/p);
  ii = length(o.arphi)*p + ii;
  if size(X,2)>0
    out.inds.X = ii:ii+size(X,2)*p-1;
    out.inds.X = reshape(out.inds.X,p,numel(out.inds.X)/p);
  else
    out.inds.X = [];
  end

  
end

function y=cfun2(d,phi)
%y = exp(-abs(d)/phi);
y = exp(-0.5*d.^2./phi.^2);
y(isnan(y))=1;

