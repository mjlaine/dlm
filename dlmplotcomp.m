function out = dlmplotcomp(dlm, t, ys, yind)
%DLMPLOTCOMP  Default plot for some DLM model components
% dlmplotcomp(dlm, t, ys, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% ys - optional y axis scale factor for plots
% yind - which series to plot, default = 1

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<2
  t = (1:size(dlm.y,1))'; % yscale for plots
end
if nargin<3
  ys = 1; % yscale for plots
end
if nargin<4
  yind = 1; % which obs column to plot
end

y = dlm.y; % obs
s = dlm.V; % obs std

p = size(dlm.F,1); % data columns
m = size(dlm.G,1); % states
mm = m/p; % states / obs columns

dlmconf = @(dlm,ind,t,ys) confband(t,ys.*dlm.x(ind,:)',ys.*dlm.xstd(:,ind));

% level and trend
ind = 1 + mm*(yind-1);
subplot(2,1,1);
dlmconf(dlm,ind,t,ys);
grid
title('level');
xlim([min(t),max(t)]);

ind = 2 + mm*(yind-1);
subplot(2,1,2);
dlmconf(dlm,ind,t,ys);
grid
title('trend');
xlim([min(t),max(t)]);

if nargout>1
  out = [];
end
