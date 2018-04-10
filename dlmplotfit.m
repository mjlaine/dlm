function out = dlmplotfit(dlm, t, ys, yind, doeb, statbox, smo)
%DLMPLOTFIT  Default plot for DLM model fit
% dlmplotfit(dlm, t, ys, yind)
% dlm - output from dlmsmo
% t - time axis for plots
% ys - optional y axis scale factor for plots
% yind - which series to plot, default = 1

% plots data with the level component and fitted values 
% with all components except possible AR

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2013/07/12 12:00:00 $

if nargin<2 | isempty(t)
  if isfield(dlm,'time')
    t = dlm.time;
  else
    t = (1:size(dlm.y,1))'; % x axis
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


% some extra inputs
if nargin<5
  doeb = 0; % do errorbar
end
if nargin<6
  statbox = 0; % add "statbox" to figures
end
if nargin<7
  smo = -1; % smoother or filter "yhat", not used yet
end

y = dlm.y; % obs
if isfield(dlm,'s')
  s = dlm.s; % obs std
else
  s = dlm.V; % obs std
end

p = size(dlm.F,1); % data columns
m = size(dlm.G,1); % states
mm = m/p; % states / obs columns

dlmconf = @(dlm,ind,t,ys) confband(t,ys.*dlm.x(ind,:)',ys.*dlm.xstd(:,ind));

% model output, excluding AR part, as the yfit will too "good"
if isfield(dlm,'options') && isfield(dlm.options,'arphi')
  nar = length(dlm.options.arphi);
else
  nar = 0;
end
if nar > 0 & 0 % now yhat is filter prediction observations
  % try to remove effect of AR
  m = size(dlm.G,1); p = size(dlm.XX,2);
  ii = 1:m-nar-p-1;
  %  ii = [ii,m-p:m];
  ii = m-nar-p:m-p;
%  yfit = ys.*(dlm.F(yind,ii)*dlm.x(ii,:))';
  yfit = ys.*(dlm.yhat(:,yind)'-dlm.F(yind,ii)*dlm.x(ii,:))';
else
  yfit = ys.*dlm.yhat(:,yind);
end

yfit = ys.*dlm.yhat(:,yind);


% data, fit,  and level
%ind = 1 + mm*(yind-1);
ind = yind; % assume cron(G,I)
hc = dlmconf(dlm,ind,t,ys);
hold on
if doeb 
  h = errorbar(t,ys.*y(:,yind),2*ys.*s(:,yind),'ok');
  if smo>=0
    hold on;plot(t,yfit,'-','color','red');hold off
  end
else
%  h = plot(t,ys.*y(:,yind),'ok',t,yfit,'-');
  h = plot(t,ys.*y(:,yind),'ok');
  if smo>=0
    hold on;plot(t,yfit,'-','color','red');hold off
  end
end
colo = [0 0.5 1];
set(h(1),'color',colo,'markerfacecolor',colo,'markersize',5);
hold off
if statbox
  % some statistics in a box
  textbox({sprintf('RMSE=%.3g',sqrt(dlm.mse(yind))),
           sprintf('MAPE=%.3g',dlm.mape(yind))});
end
  
xlim([min(t),max(t)]); % adjust xlims

% try to guess if time is matlab datetime
if t(1) > datenum(1960,1,1)
  datetick('x');
end

if isfield(dlm,'label')
  title(dlm.label);
end

if nargout>0
  out = [hc;h];
end
