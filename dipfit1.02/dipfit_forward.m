% dipfit_forward() - compute potential of a multiple dipole model
%
% Usage: 
%  >> pot = dipfit_forward(dip, elc, vol, varargin)
%
% Inputs:
%   model	dipole model for single component, can include multiple dipoles
%   elc		electrode positions
%   vol		volume conductor model
%
% Output:
%   pot		potential distribution of this model
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

% $Log: dipfit_forward.m,v $
% Revision 1.2  2003/03/03 16:51:31  roberto
% modified for posxyz/momxyz instead of dip.pos/dip.mom
% removed active and select optional inputs
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pot = dipfit_forward(model, elc, vol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a structure, which is easier to handle
if nargin>3
  optarg = struct(varargin{:});
else
  optarg = [];
end

if ~isfield(model, 'posxyz')
  error('no dipole position(s) specified');
end

if ~isfield(model, 'momxyz')
  error('no dipole moment(s) specified');
end

if isempty(model.active)
  % no active dipoles, the potential is zero
  pot = zeros(size(elc,1),1);
  return
end

% determine the total number of dipoles in the model
ndip = size(model.posxyz, 1);

% construct a reduced dipole model with only the active dipoles
dip.pos = model.posxyz(model.active,:);
dip.mom = model.momxyz(model.active,:)';

% compute the potential distribution of all active dipoles
lf  = avgref(eeg_leadfield(dip.pos, elc, vol));
pot = lf * dip.mom(:);

