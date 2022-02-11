% dipfit_relvar() - compute relative residual variance between model and data
%                   for a single component. Only active dipoles are used to 
%                   compute the model potential.
%
% Usage: 
%  >> rv = dipfit_relvar(model, data, elc, vol, varargin)
%
% Inputs:
%   data	single ICA component
%   model	dipole model for this component, can include multiple dipoles
%   elc		electrode positions
%   vol		volume conductor model
%
% Output:
%   rv		relative residual variance of dipole model
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

% $Log: dipfit_relvar.m,v $
% Revision 1.2  2003/03/03 16:52:21  roberto
% modified for posxyz/momxyz instead of dip.pos/dip.mom
%
% Revision 1.1  2003/02/24 10:05:26  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rv = dipfit_relvar(model, data, elc, vol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a structure, which is easier to handle
if nargin>4
  optarg = struct(varargin{:});
else
  optarg = [];
end

% ensure that the data is average referenced
data = avgref(data);

% compute the potential distribution of all active, unselected dipoles
pot_model = dipfit_forward(model, elc, vol);

% compute relative residual variance
rv = relvar(data, pot_model);

