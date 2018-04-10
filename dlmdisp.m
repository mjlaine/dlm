function out = dlmdisp(dlm)
%DLMDISP  print information about the DLM fit

% Marko Laine <marko.laine@fmi.fi>
% $Revision: 0.0 $  $Date: 2015/06/03 12:00:00 $


if not(isstruct(dlm)) || not(strcmp(dlm.class,'dlmfit')||strcmp(dlm.class,'dlmsmo'))
  error('works only for dlm output structure');
end

fprintf('DLM model output\n');

if strcmp(dlm.class,'dlmfit')

  fprintf('Model options\n')
  fprintf(' order: %d\n',dlm.options.order)
  if dlm.options.fullseas
    fprintf('  fullseas\n')
  else
    fprintf('  trig: %d\n',dlm.options.trig)
  end

  if not(isempty(dlm.options.arphi))
    fprintf('  AR(%d): %d\n',length(dlm.options.arphi))
  end

  if isfield(dlm,'chain')
    fprintf('MCMC: npar %d, nsimu: %d, rejected %0.3g%%\n', ...
            size(dlm.chain,2), size(dlm.chain,1),dlm.res.rejected*100);
    
  end

  fprintf('\n');

end

%fprintf('Observations %d\n',dlm.nobs);
nprint('Observations: ', '%d ',dlm.nobs);
nprint('RMSE:         ','%.3g ',sqrt(dlm.mse));
nprint('MAPE:         ','%.3g ',dlm.mape);
nprint('sigma:        ', '%.4g ',sqrt(dlm.s2));
nprint('likelihood:   ','%g ',dlm.lik);
fprintf('\n');

if nargout>0
  out=dlm;
end
function nprint(s,f,x)
fprintf('%s',s);
fprintf(f,x);
fprintf('\n');
