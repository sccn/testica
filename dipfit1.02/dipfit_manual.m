% dipfit_manual() - perform nonlinear dipole fit on one of the components
%                   to improve the initial dipole model. Only selected dipoles
%                   will be fitted. The potential of all dipoles that are active
%                   but are not selected will be subtracted from the data prior
%                   to fitting.
%
% Usage: 
%  >> dipout = dipfit_manual(model, data, elc, vol, varargin)
%
% Inputs:
%   data	single ICA component
%   model	dipole model for this component, can include multiple dipoles
%   elc		electrode positions
%   vol		volume conductor model
%
% Optional inputs are specified in key/value pairs and can be:
%   'method'		either one of 'position' / 'moment' / 'strength'
%   'constraint'	symmetry constraint structure
%
% Output:
%   dipout	result after fitting dipole model to this ICA component
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl

% Copyright (C) 2003 Robert Oostenveld, SMI/FCDC roberto@miba.auc.dk
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: dipfit_manual.m,v $
% Revision 1.5  2003/09/12 08:42:52  roberto
% changed symmetry constraint into optional input argument
%
% Revision 1.2  2003/03/03 16:51:55  roberto
% modified for posxyz/momxyz instead of dip.pos/dip.mom
%  removed active and select optional inputs
%
% Revision 1.1  2003/02/24 10:05:01  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dipout = dipfit_manual(model, data, elc, vol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a structure, which is easier to handle
if nargin>5
  optarg = struct(varargin{:});
else
  optarg = [];
end

% test the validity of the input arguments
if size(data,2)~=1
  error('data should consist of a single component only');
end
nchan = size(data,1);

if ~isfield(model, 'posxyz')
  error('no initial guess present for dipole position(s)');
end

if ~isfield(model, 'momxyz')
  error('no initial guess present for dipole moment(s)');
end

% determine the total number of dipoles in the model
ndip = size(model.posxyz, 1);

% ensure that the data is average referenced
data = avgref(data);

% process the optional arguments
if isfield(optarg, 'method')
  method = optarg.method;
else
  method = 'position';
end
if isfield(optarg, 'constraint')
  constr = optarg.constraint;
else
  constr = [];
end

if any(~ismember(model.select, model.active))
  error('all selected dipoles should be active')
end

% compute the potential distribution of all active, unselected dipoles
previous     = model.active;
model.active = setdiff(model.active, model.select);
pot_model    = dipfit_forward(model, elc, vol);	% computes potential of all active dipoles
model.active = previous;			% return to previous active dipoles

% subtract the potential distribution of all active, unselected dipoles
% the remaining data will be fitted by the selected dipoles
pot = data - pot_model;

% construct a reduced dipole model with only the selected dipoles
dip.pos = model.posxyz(model.select,:);
dip.mom = model.momxyz(model.select,:)';
nsel = length(model.select);

switch lower(method)

case 'position'
  % do nonlinear dipole fit, using the reduced input model as initial guess
  dip = eeg_dipole_fit(dip, pot, elc, vol, constr);
  % reformat the dipole moment
  dip.mom = reshape(dip.mom, 3, nsel);

case 'moment'
  % do linear estimation of rotating dipole moment
  lf = avgref(eeg_leadfield(dip.pos, elc, vol));
  dip.mom = pinv(lf) * pot;
  % reformat the dipole moment
  dip.mom = reshape(dip.mom, 3, nsel);

case 'strength'
  for i=1:nsel
    lf(:,i) = avgref(eeg_leadfield(dip.pos(i,:), elc, vol)) * dip.mom(:,i);
  end
  % estimate strength only, with fixed orientation
  strength = pinv(lf) * pot;
  for i=1:nsel
    dip.mom(:,i) = strength(i) * dip.mom(:,i);
  end

otherwise
  error('unknown method')
end

% assign fitted position and moment to the output dipole model
dipout = model;
dipout.posxyz(model.select,:) = dip.pos;
dipout.momxyz(model.select,:) = dip.mom';

% compute residual variance of complete model
dipout.rv = dipfit_relvar(dipout, data, elc, vol);

