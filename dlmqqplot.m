function out=dlmqqplot(dlm,ind)
%DLMQQPLOT normal probability plot for DLM fit residuals 

if nargin<2, ind = 1; end
% p = size(dlm.F,1); % number of series

if isfield(dlm,'resid2')
  y = dlm.resid2(:,ind);
else
  y = dlm.resid(:,ind);
end
igood = not(isnan(y));
yy = y(igood);
h=qqplot(yy);
%set(h,'markerfacecolor','black');
if isfield(dlm,'resid2')
  title('Normal probability plot for the residuals')
else
  title('Normal probability plot for the raw residuals')
end
grid on;

if nargout>0
  out=h;
end
