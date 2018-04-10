function y=acfnan(x,lagmax)
%ACFNAN Autocorrelation function
% ACF(X,maxlag)
% default maxlag is floor(10*log10(length(x)))

x = x(:)-meannan(x);
n = length(x);
if nargin<2
  lagmax = floor(10*log10(n));
  lagmax = min(lagmax, n-1);
end

y = zeros(lagmax,1);

for i=1:lagmax
  y(i) = sumnan(x(1:end-i+1).*x(i:end));
end

y = y./y(1);
