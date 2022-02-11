% dipfit_batch() - do initial batch-like dipole scan and fit to all 
%                  data components and return a dipole model with a 
%                  single dipole for each component
% 
% Usage: 
%  >> model = dipfit_batch(data, elc, vol, varargin)
%
% Inputs:
%   data    - matrix with all ICA components
%   elc     - electrode positions
%   vol     - volume conductor model
%
% Optional inputs:
%   'xgrid'   - vector with floats, grid positions along x-axis
%   'ygrid'   - vector with floats, grid positions along y-axis
%   'zgrid'   - vector with floats, grid positions along z-axis
%   'ngrid'   - integer, default number of grid positions in each direction
%   'waitbar' - ['gui'|'none'] default 'gui'.
%   'save'    - [string] save forward solutions in file given as input
%   'load'    - [string] load forward solutions in file given as input
%               this assumes that the channel location and the current 
%               volume conductor model  are identical to those used to
%               generate the file.
%
% Output:
%   model	- struct array with a dipole model for each component
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003, load/save by
%         Arnaud Delorme

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

% $Log: dipfit_batch.m,v $
% Revision 1.16  2003/12/04 18:09:05  arno
% matlab 5.3 compatibility
%
% Revision 1.15  2003/10/31 08:18:25  roberto
% changed the waitbar color into eeglab blue
%
% Revision 1.14  2003/10/14 16:08:04  roberto
% renamed grid.posxyz into pos for backward compatibility
%
% Revision 1.13  2003/10/14 15:39:47  roberto
% removed dependency on private eeg_dipole_scan
% reimplemented scanning for RV and single dipole
% reimplemented load/save of precomputed grid
%
% Revision 1.12  2003/10/09 02:33:44  arno
% same
%
% Revision 1.11  2003/10/09 02:30:14  arno
% same
%
% Revision 1.10  2003/10/09 02:29:31  arno
% save
%
% Revision 1.9  2003/10/09 02:28:08  arno
% partial save abord
%
% Revision 1.8  2003/10/09 02:18:18  arno
% cell array input problem
%
% Revision 1.7  2003/10/09 01:19:22  arno
% implementing additional options
%
% Revision 1.6  2003/06/30 01:42:18  arno
% waitbar -> feedback
%
% Revision 1.5  2003/06/30 01:39:58  arno
% waitbar argument
%
% Revision 1.4  2003/03/06 15:57:37  roberto
% *** empty log message ***
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = dipfit_batch(data, elc, vol, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert the optional arguments into a structure, which is easier to handle
if nargin>3
    for index = 1:length(varargin)
        if iscell(varargin{index})
            varargin{index} = {varargin{index}};
        end;
    end;
    optarg = struct(varargin{:});
else
  optarg = [];
end

% process the optional arguments
if isfield(optarg, 'ngrid')
  ngrid = optarg.ngrid;
else
  ngrid = 15;
end
if isfield(optarg, 'xgrid')
  xgrid = optarg.xgrid;
else
  xgrid = [];	% will be assigned the default value later on
end
if isfield(optarg, 'ygrid')
  ygrid = optarg.ygrid;
else
  ygrid = [];	% will be assigned the default value later on
end
if isfield(optarg, 'zgrid')
  zgrid = optarg.zgrid;
else
  zgrid = [];	% will be assigned the default value later on
end
if isfield(optarg, 'waitbar') & strcmp(optarg.waitbar, 'none')
    waitbarflag = 0;
else
    waitbarflag = 1;
end;
loadflag = isfield(optarg, 'load');
saveflag = isfield(optarg, 'save');

if ~loadflag
  % create a grid
  if isempty(xgrid) | isempty(ygrid) | isempty(zgrid) 
    % create defaults for the grid along each axis
    if isfield(vol, 'r')
      % spherical model, use innermost sphere diameter to define scan space
        brain_radius = min(vol.r);
        if isfield(vol, 'o')
          % the center of the sphere is not in the origin
          brain_center = vol.o;
        else
          brain_center = [0 0 0];
        end
      xmin = brain_center(1) - brain_radius; xmax = brain_center(1) + brain_radius;
      ymin = brain_center(2) - brain_radius; ymax = brain_center(2) + brain_radius;
      zmin = brain_center(3) - brain_radius; zmax = brain_center(3) + brain_radius;
      clear brain_radius brain_center
    elseif isfield(vol, 'bnd')
      % BEM model, use the dimensions of the brain surface to define scan space
      % FIXME, implement this at a later time
    end
    if isempty(xgrid), xgrid = linspace(xmin, xmax, ngrid); end
    if isempty(ygrid), ygrid = linspace(ymin, ymax, ngrid); end
    if isempty(zgrid), zgrid = linspace(zmin, zmax, ngrid); end
    clear xmin xmax ymin ymax zmin zmax
  end
  % create the regular grid with dipole positions
  [X, Y, Z] = meshgrid(xgrid, ygrid, zgrid);
  grid.pos = [X(:), Y(:), Z(:)];
  % select only those grid positions that are inside the brain
  [inside, outside] = find_inside_vol(grid.pos, vol);
  grid.pos = grid.pos(inside,:);
else
  % read the grid and precomputed inverted leadfield from an external file
  load(optarg.load);
end

% assign the output, start with an empty model for each component
for comp=1:size(data,2)
  model(comp).posxyz = [];
  model(comp).momxyz = [];
  model(comp).rv     = inf;
end

if waitbarflag
  % get the defaults for the GUI
  icadefs;
  try 
      waitbarhandle = waitbar(0, 'Scanning ...', 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 'Color', BACKEEGLABCOLOR);
  catch % for Matlab 5.3
      waitbarhandle = waitbar(0, 'Scanning ...');
  end;   
end

% compute the residual variance for each component and grid location
nchan  = size(data,1);
ndip   = size(grid.pos,1);
data   = avgref(data);
sumsqr = sum(data.^2,1);
for dip=1:size(grid.pos,1)
  if waitbarflag
    waitbar(dip/ndip, waitbarhandle);
  end

  if loadflag
    % use the precomputed inverse leadfield that was stored along with the grid
    lfi = grid.lfi{dip};
  else
    lf  = avgref(eeg_leadfield(grid.pos(dip,:), elc, vol));
    lfi = pinv(lf);
  end

  if saveflag
    % store the computed inverse leadfield together with the grid so that it can be saved
    grid.lfi{dip} = lfi;
  end

  rv = sum(((eye(nchan) - lf*lfi) * data) .^2, 1) ./ sumsqr;
  for comp=1:size(data,2)
    if rv(comp) < model(comp).rv 
      % this position has a lower residual variance that the present model
      model(comp).rv     = rv(comp);
      model(comp).posxyz = grid.pos(dip,:);
      model(comp).momxyz = lfi*data(:,comp);	% recompute dipole moment
      model(comp).momxyz = model(comp).momxyz';	% and make it a row vector
    end
  end
end

if saveflag
  % save the grid together with the precomputed inverse leadfields to an external file
  save(optarg.save, 'grid');
end

if waitbarflag & ishandle(waitbarhandle)
  delete(waitbarhandle);
end

