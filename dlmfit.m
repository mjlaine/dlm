function out = dlmfit(y,s,wdiag,x0,C0, X, options)
%DLMFIT Fit DLM time series model 
% Fits dlm time series model with local level, trend, seasonal, and proxies
% out = dlmfit(y,s,wdiag,x0,C0, X, options)
% Input:
% y time series, n*p
% s obs uncertainty, n*p or 1*1
% w sqrt of first diagonal entries of the model error matrix W
% x0 initial state m*1 (optional, i.e. can be empty)
% C0 initial state uncertainty covariance matrix m*m (optional)
% X covariate variables n*q (optional)
% options structure
%   order
%   fullseas
%   trig
%   ns
%   opt
%   mcmc
% Output:
% out.

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin < 6
  X=[]; % covariates, aka proxies
end

if nargin < 7
  options.trig = 0;
  options.opt  = 0;
  options.mcmc = 0;
end

% compatibility for older version
if isfield(options,'fullseas'), options.seas=options.fullseas; end

if not(isfield(options,'order')), options.order=1; end
if not(isfield(options,'trig')),  options.trig=0; end
if not(isfield(options,'seas')),  options.seas=0; end
if options.trig>0, options.seas=0; end
if not(isfield(options,'fullseas')), options.fullseas=options.seas; end

if not(isfield(options,'ns')), options.ns=12; end % number of seasons

if not(isfield(options,'mcmc')), options.mcmc=0; end
if not(isfield(options,'opt')),  options.opt=0; end
if not(isfield(options,'maxfuneval')),  options.maxfuneval=400; end

% fit V factor?
if not(isfield(options,'fitv')),  options.fitv=0; end % see also options.vcv

if not(isfield(options,'winds')),  options.winds=[]; end % 
if not(isfield(options,'logscale')),  options.logscale=1; end % 

if not(isfield(options,'spline')),  options.spline=0; end % test

% These are the defaults used in the dlm acp paper
if not(isfield(options,'varcv')), options.varcv=[1 1 1 1]; end
if not(isfield(options,'vcv')), options.vcv=0.5; end % V factor prior CV
if not(isfield(options,'nsimu')), options.nsimu=5000; end

[G,F] = dlmgensys(options); % generate system matrises

[p,m] = size(F); % series states
n = length(y); % nobs

% add covariates to system matrix G
if not(isempty(X))
  kk = size(X,2);
  G(m+1:m+kk,m+1:m+kk) = eye(kk);
  m = m+kk;
end

V  = ones(n,p).*s;    % V is matrix of std's
W  = zeros(m,m);      % model error

% input wdiag has diagonal std of W
for i=1:length(wdiag); W(i,i) = wdiag(i).^2; end

% spline
if options.order==1 & options.spline
  W(1:2,1:2) = wdiag(2).^2 .* [1/3 1/2;1/2 1];
end


% try to find sensible initial values (FIX THIS)
if nargin < 4 || isempty(x0)
  x0 = zeros(m,1);   % initial
  x0(1) = meannan(y(1:ceil(options.ns),:)); % how about p>1 (not yet) assume x0(1) is the level
%  x0(2) = (y(13)-y(1))/12; % crude estimate
%  if options.trig==0
%    x0(3) = sumnan(-y(1:11,:)); % -,,-
%  end
end
if nargin<5||isempty(C0)
  C0diag = ones(1,m)*(abs(x0(1))*0.5).^2;
  C0diag(C0diag==0) = 1e+7;
  C0 = diag(C0diag);
end
% fit, and refit to iterate the initial values
out = dlmsmo(y,F,V,x0,G,W,C0, X, 0);
x0 = out.x(:,1);
C0 = 100*squeeze(out.C(:,:,1));
out = dlmsmo(y,F,V,x0,G,W,C0, X);

%% optimization and MCMC calculations

if options.opt
  if isfield(options,'fitfun')
    % use user given fitfun that returs dlm object
    out = options.fitfun(out,options);    
  else
    % default optimization for a V and W
    out = dlm_dooptim(out,options);
  end
end

if options.mcmc
  if isfield(options,'mcmcfun')
    out = options.mcmcfun(out,options);
  else
    % default MCMC for W and V
    out = dlm_domcmc(out,options);
  end
end

out.s = s; % save original obs std also
out.options = options;
out.class = 'dlmfit';

%% some helper functions

function out = dlm_dooptim(dlm,options)
% this optimizes some parameters defining the DLM model
%  the costfun is at the end

if options.fitv > 0; 
  vinds = 1;
else
  vinds = [];
end

% winds = map from diag(W) to optimized parameter
winds = options.winds;
nw = max(winds);

if isfield(options,'fitar') && options.fitar>0
  ginds = findarinds(options);
else
  ginds = [];
end

zfun = @(x) dlm_costfun(x,dlm,vinds,winds,ginds); 

W = dlm.W;
V = dlm.V;
G = dlm.G;

if length(vinds)>0;
  v0 = 1;
else 
  v0 = [];
end 

w0 = zeros(nw,1);
for i=1:length(w0)
  ii=find(winds==i,1);
  w0(i)=sqrt(W(ii,ii));
end

if length(ginds)>0
  g0 = options.arphi(:);
else
  g0 = [];
end

w00 = [log([v0;w0]);g0];
inds = [ones(size(v0));ones(size(w0))*2;ones(size(g0))*3];

oo = optimset('disp','iter','maxfuneval',options.maxfuneval);
wopt = fminsearch(zfun,w00,oo);

woptv = exp(wopt(inds==1))
woptw = exp(wopt(inds==2))
woptg = wopt(inds==3)

if length(vinds)>0;
  V = V.*woptv;
end

for i=1:length(winds);
  if winds(i)>0;
    W(i,i) = woptw(winds(i)).^2;
  end;
end

arinds = findarinds(options);
if length(ginds)>0
  G(arinds,arinds(1)) = woptg;
end

out = dlmsmo(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX);

function out = dlm_costfun(x,dlm,vinds,winds,ginds)
% Helper function for dlm parameter fitting
W = dlm.W;
V = dlm.V;
G = dlm.G;

nv = length(vinds);
nw = max(winds);
ng = length(ginds);

if nv>0;
  V=V.*exp(x(1));
end
for i=1:length(winds)
  if winds(i)>0;
    W(i,i)=exp(x(winds(i)+nv)).^2;
  end
end
if ng>0
  G(ginds,ginds(1)) = x(nv+nw+1:end);
end

if exist('dlmmex') == 3
  % mex boosted likelihood calculations
  out = dlmmex(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX);
else
  out = getfield(dlmsmo(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX,0,0),'lik');
end

function out = dlm_domcmc(dlm,options)
% this does MCMC for model parameters

% fit a coefficient for obs error
if options.fitv > 0; 
  vinds = 1;
else
  vinds = [];
end

% winds = map from diag(W) to optimized parameter
winds = options.winds;
nw = max(winds);

if isfield(options,'fitar') && options.fitar>0
  ginds = findarinds(options);
else
  ginds = [];
end

% save current values
W = dlm.W;
V = dlm.V;
G = dlm.G;

% initial values
if length(vinds)>0;
  v0 = 1;
else 
  v0 = [];
end 

w0 = zeros(nw,1);
for i=1:length(w0)
  ii=find(winds==i,1);
  w0(i)=sqrt(W(ii,ii));
end
w0ini = w0;
% use option for initial values, must be of length w0
if isfield(options,'varini')
  if length(options.varini) == length(w0)
    w0ini = options.varini(:);
  else
    error('options.varini length does not match')
  end
end

if length(ginds)>0
  g0 = options.arphi(:);
else
  g0 = [];
end

w00 = [log([v0;w0]);g0];
inds = [ones(size(v0));ones(size(w0))*2;ones(size(g0))*3];

if options.logscale
  ffun=@(x)log(x);
  bfun=@(x)exp(x);
  parmin = -Inf;
else
  ffun=@(x)(x);
  bfun=@(x)(x);
  parmin = 0;
end

npar = length(w00);
p = cell(0);
i1 = 0;
for i=1:length(v0)
  i1 = i1+1;
  p{i1} = {'Vfact',ffun(v0(i)),parmin,Inf,NaN,options.vcv(i)};
end
for i=1:length(w0)
  i1 = i1+1;
  ii = find(winds==i,1);
%  p{i1} = {sprintf('w%d',ii),ffun(w0(i)),parmin,Inf,NaN,options.varcv(i)};
  p{i1} = {sprintf('w%d',ii),ffun(w0ini(i)),parmin,Inf,ffun(w0(i)),options.varcv(i)};
end
for i=1:length(g0)
  i1 = i1+1;
  % bound [0,1] for AR rho is fixed here for now
  p{i1} = {sprintf('g%d',i),g0(i),0,1,NaN,options.gcv(i)};
end

m = struct; o = struct;
if options.logscale
  m.ssfun = @(th,d) fitfun_mcmc(th,dlm,vinds,winds,ginds);
  m.priorfun = @(th,m,s) sum(((th-m)./s).^2); % normal prior
  m.priortype = 1; % 
else
  m.ssfun = @(th,d) error('not defined yet');
  m.priorfun = @(th,m,s) sum((log(th./m)./s).^2); % lognormal prior
  m.priortype=-1; % for lognormal prior in densplot
end
o.qcov = eye(length(w00));
i1 = 0;
for i=1:length(v0)
  i1 = i1+1;
  o.qcov(i1,i1) = 0.1; % proposal for V fact fixed here
end
for i=1:length(w0)
  i1 = i1+1;  
  o.qcov(i1,i1) = abs(p{i1}{2})*0.002; %  W diag paramrs proposal
  if o.qcov(i1,i1) == 0, o.qcov(i1,i1) = 1; end
end
for i=1:length(g0)
  i1 = i1+1;  
  o.qcov(i1,i1) = 0.01; % propsal for AR rhos
end

o.nsimu = options.nsimu;
o.method = 'am';
o.adaptint = 100;
o.initqcovn = 200;
o.verbosity = 1;
[res,chain,~,sschain] = mcmcrun(m,[],p,o);

wopt = mean(chain(fix(size(chain,1)/2):end,:)); % use the last half

woptv = bfun(wopt(inds==1))
woptw = bfun(wopt(inds==2))
woptg = wopt(inds==3)

if length(vinds)>0;
  V = V.*woptv;
end

for i=1:length(winds);
  if winds(i)>0;
    W(i,i) = woptw(winds(i)).^2;
  end;
end

arinds = findarinds(options);
if length(ginds)>0
  G(arinds,arinds(1)) = woptg;
end

out = dlmsmo(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX);

if options.logscale % scale chain back to stds
  ii = (inds ==1|inds==2);
  chain(:,ii) = bfun(chain(:,ii));
  res.prior(ii,1)=exp(res.prior(ii,1));
  res.priorfun = @(th,m,s) sum((log(th./m)./s).^2);
  res.priortype=-1;
end

out.res = res;
out.chain = chain;
out.sschain = sschain;
out.winds = winds;
out.vinds = vinds;
out.ginds = ginds;


function out = fitfun_mcmc(x,dlm,vinds,winds,ginds)
% Helper function for dlm parameter fitting
W = dlm.W;
V = dlm.V;
G = dlm.G;

nv = length(vinds);
nw = max(winds);
ng = length(ginds);


if nv>0;
  V=V.*exp(x(1));
end
for i=1:length(winds)
  if winds(i)>0;
    W(i,i)=exp(x(winds(i)+nv)).^2;
  end
end
if ng>0
  G(ginds,ginds(1)) = x(nv+nw+1:end);
end

if exist('dlmmex') == 3
  % mex boosted likelihood calculations
  out = dlmmex(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX);
else
  out = getfield(dlmsmo(dlm.y,dlm.F,V,dlm.x0,G,W,dlm.C0,dlm.XX,0,0),'lik');
end

function i = findarinds(options)
% find AR indeces
if not(isfield(options,'arphi'));
  i = [];
  return;
end
nar = length(options.arphi);
i = 1;
i = i + options.order + 1;
if options.fullseas==1
  i = i + 11;
elseif options.trig>0
  i = i+min(options.ns-1,options.trig*2);
end

i = i:i+nar-1;
