function out = dlmplotdiag(dlm, t, ys, yind)
%DLMPLOTDIAG  Default plot for DLM model diagnostics
% dlmplotdiag(dlm, t, ys, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% ys - optional y axis scale factor for plots (not used)
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

% acf and qqplot
subplot(2,1,1)
dlmacfplot(dlm,yind);
subplot(2,1,2)
dlmqqplot(dlm,yind);

if nargout>1
  out = [];
end
