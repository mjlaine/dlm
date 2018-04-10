function out = dlmplot(dlm, t, ys, yind)
%DLMPLOT  Default plot for DLM model fit
% dlmplot(dlm, t, ys, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% ys - optional y axis scale factor for plots
% yind - which series to plot, default = 1

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<2 | isempty(t)
  if isfield(dlm,'time') 
    t = dlm.time;
  else
    t = (1:size(dlm.y,1))';
  end
end
if nargin<3 | isempty(ys)
  if isfield(dlm,'ys') 
    ys = dlm.ys;
  else
    ys = 1; % yscale for plots
  end
end
if nargin<4
  yind = 1; % which obs column to plot
end

doeb = 0; % do errorbar
fex = 10*(yind-1); % figure offset

% data and fit
figure(fex+1); clf
dlmplotfit(dlm, t, ys, yind, doeb, 0, 0);
% level and trend
figure(fex+2); clf
dlmplotcomp(dlm, t, ys, yind);
% acf and qqplot
figure(fex+3); clf
dlmplotdiag(dlm, t, ys, yind);

if nargout>1
  out = [];
end
