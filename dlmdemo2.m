%% DLM demo no 2, testing DLM functions
% Generate seasonal data and fit DLM model with 12 seasons.

%% 
% Data
t = (1:1:(6*12))'; % time 
n = length(t);
s = 0.1; % obs error std 
y = 2 + 0.5*(t-0)/12.*((t-0)<2*12)    + ...
        0.5*(2*12)/12.*((t-0)>=2*12)  + ...
        0.8*sin((t-0)/12*2*pi) + ...
     randn(n,1)*s;

%%
figure(1); clf
plot(t,y,'o-'); xlabel('time'); ylabel('y'); title('observations');grid;

%% 
% Model error diagonal, first non zero values.
w123 = [0.0000,3.1e-5,1e-8,1e-8,1e-8,1e-8];

%%
% Options for |dlmfit| function. Use harmonic functions with 2
% components and 12 seasons. Optimize some variance parameters using maximum likelihood.
clear options
options.ns = 12;
options.trig = 2;
options.opt = 1;
options.winds = [0 1 2 2 2 2];

%% 
% DLM fit.
out = dlmfit(y,s,w123,[],[],[],options);
%%
% Default plots.
dlmplot(out);

