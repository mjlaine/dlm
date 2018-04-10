function out = dlmplotresid(dlm, t, yind)
%DLMPLOTRESID  Plot for DLM residuals
% dlmplotfit(dlm, t, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% yind - which series to plot, default = 1

% plots data with the level component and fitted values 
% with all components except possible AR

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<2
  if isfield(dlm,'time')
    t = dlm.time;
  else
    t = (1:size(dlm.y,1))'; 
  end
end
if nargin<3
  yind = 1; % which obs column to plot
end

colo = [0 0.5 1];
siz = 4;
setplot = @(h) set(h,'color',colo,'markersize',siz,'markerfacecolor',colo);

subplot(3,1,1)
h1=plot(t,dlm.resid0(:,yind),'o-');setplot(h1);grid on
title('raw smoother residuals')
xlim([t(1),t(end)]);

subplot(3,1,2)
h2=plot(t,dlm.resid(:,yind),'o-');setplot(h2);grid on
title('scaled smoother residuals')
xlim([t(1),t(end)]);

subplot(3,1,3)
h3=plot(t,dlm.resid2(:,yind),'o-');setplot(h3);grid on
title('scaled prediction residuals')
xlim([t(1),t(end)]);

if nargout>0
  out=[h1;h2;h3];
end
