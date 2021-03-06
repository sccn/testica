% pop_dipoles() - export components for localisation in BESA
%
% Usage:
%   >> pop_dipoles( maparray ); % gui pops-up
%   >> pop_dipoles( maparray, components, 'key', 'val', ...);
%   >> pop_dipoles( EEG ); % gui pops-up
%   >> pop_dipoles( EEG, components, 'key', 'val', ...);
%
% Inputs:
%   maparray   - Input scalp map array of size (electrodes, components)
%   components - [integer vector] Component indices to export
%   EEG        - An eeglab() dataset (with associated location file)
%                whose components will be dipole-localized.
%
% Optional inputs:
%   'dipoles'  - [integer|integer array] Number of dipoles per component. 
%                Either a single value or a value array (1 per component).
%                {Default 1}.
%   'elpfile'  - External '.elp' BESA file for electrode location. 
%   'elecfile' - All-electrode location file compatible with readlocs().
%                Electrode positions are converted into BESA spherical coord.
%   'elecshrink' - Shrink factor used for polar-coordinate transformation.
%   'averef'   - [integer integer] Average-reference localization, if any. 
%                Azimut then horizontal angle as in BESA. Ex: "0 0" for Cz.
%
% See also: besaplot(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C)
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

% $Log: pop_dipoles.m,v $
% Revision 1.1  2003/02/14 00:23:48  arno
% Initial revision
%

function com = pop_dipoles( winvarray, components, varargin )

com='';
if nargin < 1
   help pop_dipoles;
   return;
end;

if nargin < 2

	% popup window parameters
	% -----------------------
	components = [];
	promptstr    = { 'Component numbers to localize:' ...
					 strvcat('Number of dipoles per component', ...
							 '(a single value or one value per component)') ...
					 strvcat('External electrode locations file', ...
							 '(If using an eeglab() dataset, location file will be generated automatically)', ...
							 '(For other file types, call this function from the command line)') ...
					 strvcat('Optional location of the common-reference electrode', ...
							 '(azimuth then horizontal 3-D ngle as in BESA. Ex "0 0" for Cz)', ...
							 '(for average reference data, leave this field blank)'), ...
					 strvcat('Volume conductor model') };

	inistr       = { '1:10' ...
					 '1' ...
					 '' ...
					 '' ...
                     'standard.vol' };

    result       = inputdlg2( promptstr, 'Localize equivalent dipoles -- pop_dipoles()', ...
                              1,  inistr, 'pop_dipoles');
	if length(result) == 0 return; end;
	components     = eval( [ '[' result{1} ']' ] );
	g.dipoles      = eval( [ '[' result{2} ']' ] );
	g.elpfile      = result{3};
	g.averef       = result{4};
   	g.volume       = result{5};
 
    if isstruct(winvarray)
		g.elecfile = winvarray.chanlocs;
		winvarray  = winvarray.icawinv;
	else 
		g.elecfile = [];
	end;	
	g.elecshrink = 1;
    g
    
	% interface for scanning components
	% ---------------------------------
    geomhoriz = { [1 1] [1] [1 0.5 2] [1] [1 1 1 2 2 2] [1 0.7 0.2 0.7 0.2 2 2 2] [1 0.7 0.2 0.7 0.2 2 2 2] ...
                  [1 0.7 0.2 0.7 0.2 2 2 2]  [1 0.7 0.2 0.7 0.2 2 2 2] [1] [1 1 1 1]};
    geomvert  = [2 1 1 1 1 1 1 1 1 1 1];
	textgui      = { ...
                     { 'style' 'pushbutton'  'string' 'Fit dipoles for all components' } ...   
                     { 'style' 'pushbutton'  'string' 'Show all dipoles (3-D)' } ...   
                     { } ...
                     { 'style' 'text'    'string' 'Select ICA component' } ...
                     { 'style' 'listbox' 'string' int2str(components') } ...
                     { 'style' 'text'    'string' 'Residual variance:' 'tag' 'resvar' } ...
                     { } ...
                     { 'style' 'text'    'string' 'dipole' } ...
                     { 'style' 'text'    'string' 'active' } ...
                     { 'style' 'text'    'string' 'select' } ...
                     { 'style' 'text'    'string' 'constraint' } ...
                     { 'style' 'text'    'string' 'position' } ...
                     { } ...
                     { 'style' 'text'        'string' '#1' } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'listbox'     'string' strvcat('none|symmetry to ?|...') } ...
                     { 'style' 'edit'        'string' ' ' } ...
                     { 'style' 'pushbutton'  'string' 'Scan position' } ...   
                     ...
                     { 'style' 'text'        'string' '#2' } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'listbox'     'string' strvcat('none|symmetry') } ...
                     { 'style' 'edit'        'string' ' ' } ...
                     { 'style' 'pushbutton'  'string' 'Scan position' } ...                   
                     ...
                     { 'style' 'text'        'string' '#3' } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'listbox'     'string' strvcat('none|symmetry') } ...
                     { 'style' 'edit'        'string' ' ' } ...
                     { 'style' 'pushbutton'  'string' 'Scan position' } ...                   
                     ...
                     { 'style' 'text'        'string' '#4' } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'checkbox'    'string' ' ' } { } ...
                     { 'style' 'listbox'     'string' strvcat('none|symmetry') } ...
                     { 'style' 'edit'        'string' ' ' } ...
                     { 'style' 'pushbutton'  'string' 'Scan position' } ...                   
                     { } ...
                     { 'style' 'pushbutton'  'string' 'Fit dipoles' } ...   
                     { 'style' 'pushbutton'  'string' 'Fit dipoles moments' } ...   
                     { 'style' 'pushbutton'  'string' 'Show dipoles (3-D)' } ...   
                     { 'style' 'pushbutton'  'string' 'Show potential map' } ...   
                   };

    result = inputgui(geomhoriz, textgui, 'pophelp(''pop_dipole'');', 'pop_dipole', [], 'normal', geomvert);
end;



