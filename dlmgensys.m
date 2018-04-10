function [G,F] = dlmgensys(options)
%DLMGENSYS generate system matrices G and F for a DLM model
% dlmgensys(options)
% options.order    - order of the local polynomial trend component (1)
% options.fullseas - to use full seasonal component or not (0)
% options.trig     - how many seasonal harmonics (0)
% options.ns       - how many seasons (default = 12)
% opitions.arphi   - AR coefficients ([])

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<1
  options = struct();
end

% check options field
if not(isfield(options,'order')),     options.order=1; end
if not(isfield(options,'fullseas')),  options.fullseas=0; end
if not(isfield(options,'trig')),      options.trig=0; end
if not(isfield(options,'ns')),        options.ns=12; end

if not(isfield(options,'arphi')),     options.arphi=[]; end

if options.trig>options.ns/2
  error('trig must be between 0 and ns/2')
end

% utility to combine matrices
stack = @(a,b)[a,zeros(size(a,1),size(b,2));zeros(size(b,1),size(a,2)),b];
% Generate harmonic component
ns = options.ns;
harm = @(k)[cosd(360/ns*k),sind(360/ns*k);-sind(360/ns*k),cosd(360/ns*k)];
% Generate polynomial trend component
trend = @(o) diag(ones(o+1,1))+diag(ones(o,1),1);

% Local polynomial trend
Gt = trend(options.order);
Ft = [1, zeros(1,options.order)];
if isempty(Gt) , Ft = []; end

% Seasonal term
if options.fullseas == 1
  Gs = diag(ones(options.ns-2,1),-1);
  Gs(1,:) = -1;
  Fs = [1,zeros(1,options.ns-2)];
else
  Gs = [];
  Fs = [];
  for i = 1:options.trig
    Gs = stack(Gs,harm(i));
    Fs = [Fs,[1 0]];
  end
  % if ns even and trig=ns/2, then remove the last element as it is redundant
  if options.trig == options.ns/2
    Gs = Gs(1:end-1,1:end-1);
    Fs = Fs(1:end-1);
  end
end

% AR terms (experimental)
arphi=options.arphi(:);
nar = length(arphi);
if nar > 0
  Gar = [arphi,[eye(nar-1);zeros(1,nar-1)]];
  Far = [1,zeros(1,nar-1)];
else
  Gar = [];
  Far = [];
end

% combine trend and seasonal
G = stack(Gt,Gs);
F = [Ft,Fs];

% add AR
G = stack(G,Gar);
F = [F,Far];

