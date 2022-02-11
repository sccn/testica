function [dip] = eeg_dipole_fit(dip, dat, elc, vol, constr)

% EEG_DIPOLE_FIT performs an equivalent current dipole fit with a 
% single dipole or a small number of dipoles.
%
% [dip] = eeg_dipole_fit(dip, dat, elc, vol, constr)

% Copyright (C) 2003, Robert Oostenveld
%
% $Log: eeg_dipole_fit.m,v $
% Revision 1.6  2003/09/12 08:45:41  roberto
% added error if maximum number of iterations is exceeded
%
% Revision 1.5  2003/09/02 13:01:18  roberto
% implemented constrained dipole fitting
%
% Revision 1.4  2003/06/16 10:03:31  roberto
% added check in error function for gui interrupt request (for eeglab)
%
% Revision 1.3  2003/03/13 13:41:00  roberto
% fixed bug in e/meg_error_func in assignment of Nchan
%
% Revision 1.2  2003/03/11 14:45:36  roberto
% updated help and copyrights
%

if nargin<5
  constr = [];
end

% ensure that the reference of the potential data is the same everywhere
dat = avgref(dat);

% reformat the position parameters in case of multiple dipoles
param = reshape(dip.pos', 1, prod(size(dip.pos)));

% reduce the number of parameters to be fitted according to the constraints
if ~isempty(constr)
  param = param(constr.reduce);
end

if exist('fminunc')
  % do non-linear optimization of the dipole parameters using Optimization toolbox fminunc()
  options = optimset('TolFun',1e-9,...
                     'TypicalX',ones(size(param)),...
                     'LargeScale','off',...
                     'MaxIter',100,...
                     'HessUpdate','bfgs',...
                     'Display','iter');
  [param, fval, exitflag, output] = fminunc(@dipfit_error, param, options, dat, elc, vol, constr);
else
  % do non-linear optimization of the dipole parameters using standard Matlab fminsearch()
  options = optimset('MaxIter',100,...
                     'Display','iter');
  [param, fval, exitflag, output] = fminsearch(@dipfit_error, param, options, dat, elc, vol, constr);
end

if exitflag==0
  error('Maximum number of iterations exceeded before reaching the minimum, I suggest that you try with another initial guess.')
end

% do linear optimization of dipole moment parameters
[err, dip.mom] = dipfit_error(param, dat, elc, vol, constr);

% expand the number of parameters according to the constraints
if ~isempty(constr)
  param = constr.mirror .* param(constr.expand);
end

% reformat the position parameters in case of multiple dipoles
dip.pos = reshape(param, 3, length(param)/3)';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIPFIT_ERROR computes the error between measured and model data
% and can be used for non-linear fitting of dipole position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, mom] = dipfit_error(param, dat, elc, vol, constr);

% flush pending graphics events, ensure that fitting is interruptible
drawnow;
if ~isempty(get(0, 'currentfigure')) & strcmp(get(gcf, 'tag'), 'stop')
  % interrupt the fitting
  close;
  error('USER ABORT');
end;

% expand the number of parameters according to the constraints
if ~isempty(constr)
  pos = constr.mirror .* param(constr.expand);
else
  pos = param;
end

Ndipoles = length(pos)/3;
if isfield(vol, 'tra')
  Nchans   = size(vol.tra,1);
else
  Nchans   = size(dat,1);
end

% this ensures that the position is a row-vector
pos = pos(:)';

if mod(length(pos),3)
  error('invalid number of positions (should be Nx3)');
end

% construct the leadfield matrix for all dipoles
lf = zeros(Nchans, 3*Ndipoles);
for i=1:Ndipoles
  sel = (3*(i-1)+1):(3*(i-1)+3);
  lf(:, sel) = eeg_leadfield(pos(sel), elc, vol);
end

% ensure that the reference of the potential data is the same everywhere
lf = avgref(lf);

% compute the optimal dipole moment and the model error
mom = pinv(lf)*dat;
dif = dat - lf*mom;
err = sum(dif(:).^2) / sum(dat(:).^2);

