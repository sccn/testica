% pop_dipfit_settings() - interactively change the global settings for dipole fitting
%
% Usage:
%   >> OUTEEG = pop_dipfit_settings ( INEEG ); % pop up window
%   >> OUTEEG = pop_dipfit_settings ( INEEG, 'key1', 'val1', 'key2', 'val2' ... )
%
% Inputs:
%   INEEG	input dataset
%
% Optional inputs:
%   'radii'        - [float array] radii values for the model above.
%                    Default is [] (uses default)
%   'conductances' - [float array] conductance values for the model above.
%                    Default is [] (uses default)
%   'electrodes'   - [integer array] indices of electrode to include
%                    in model. Default: all.
%
% Outputs:
%   OUTEEG	output dataset
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%         Arnaud Delorme, SCCN, La Jolla 2003

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

% $Log: pop_dipfit_settings.m,v $
% Revision 1.28  2005/03/17 00:34:57  arno
% fix selection
%
% Revision 1.27  2004/01/07 17:02:48  scott
% See List -> List
%
% Revision 1.26  2004/01/07 17:01:32  scott
% ... -> See List
%
% Revision 1.25  2004/01/07 17:00:21  scott
% added to Note
%
% Revision 1.24  2003/12/04 18:25:43  arno
% shell conductance unit
%
% Revision 1.23  2003/10/31 18:06:07  arno
% more edits
%
% Revision 1.22  2003/10/31 17:54:00  arno
% select channels -> omit channels
%
% Revision 1.21  2003/10/30 02:14:51  arno
% gui typo
%
% Revision 1.20  2003/10/29 22:43:37  arno
% wording
%
% Revision 1.19  2003/10/29 02:36:38  arno
% removing 2 lines of GUI
%
% Revision 1.18  2003/10/15 14:47:34  roberto
% removed urchanlocs from electrode projection part
%
% Revision 1.17  2003/10/14 15:56:38  roberto
% before projecting electrodes towards the skin, store the original in urchanlocs
%
% Revision 1.16  2003/08/08 16:57:10  arno
% removing normsphere
%
% Revision 1.15  2003/08/04 22:10:13  arno
% adding warning backtrace
%
% Revision 1.14  2003/08/04 22:03:32  arno
% adding normsphere option
%
% Revision 1.13  2003/08/01 13:50:51  roberto
% changed the gui, implemented fit-electrodes-to-sphere, modifications in vol.r/c/o are accepter
%
% Revision 1.12  2003/07/01 22:11:59  arno
% removing debug message
%
% Revision 1.11  2003/06/30 02:11:37  arno
% debug argument check
%
% Revision 1.10  2003/06/30 01:20:58  arno
% 4sphere -> 4spheres
%
% Revision 1.9  2003/06/30 01:18:46  arno
% copying new version
%
% Revision 1.5  2003/06/16 15:32:51  arno
% reprograming interface, programing history
%
% Revision 1.4  2003/03/12 10:32:50  roberto
% added 4-sphere volume model similar to BESA
%
% Revision 1.3  2003/03/06 15:58:28  roberto
% *** empty log message ***
%
% Revision 1.1  2003/02/24 10:06:08  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_settings ( EEG, varargin )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help pop_dipfit_settings;
   return;
end;

OUTEEG = EEG;
com = '';

% get the default values and filenames
dipfitdefs;

if nargin < 2
    % define the callbacks for the buttons
    cb_selectelectrodes = 'tmp = select_channel_list({EEG.chanlocs.label}, eval(get(findobj(gcbf, ''tag'', ''elec''), ''string''))); set(findobj(gcbf, ''tag'', ''elec''), ''string'',[''[''  num2str(tmp) '']''])'; % did not work
    cb_selectelectrodes = 'set(findobj(gcbf, ''tag'', ''elec''), ''string'', int2str(pop_chansel({EEG.chanlocs.labels})));';
    % cb_selectvolume     = '[fname, pname] = uigetfile(''*.*'', ''Load volume conductor''); fname = fullfile(pname, fname); tmp = strvcat(get(findobj(gcbf, ''tag'', ''vol''), ''string''), fname); set(findobj(gcbf, ''tag'', ''vol''), ''string'', tmp); set(findobj(gcbf, ''tag'', ''vol''), ''value'', size(tmp,1)); set(gcbf, ''userdata'', tmp)';
    cb_volmodel = [ 'tmpdat = get(gcbf, ''userdata'');' ... 
                    'tmpind = get(gcbo, ''value'');' ... 
                    'set(findobj(gcbf, ''tag'', ''radii''),   ''string'', num2str(tmpdat{tmpind}.r,3));' ...
                    'set(findobj(gcbf, ''tag'', ''conduct''), ''string'', num2str(tmpdat{tmpind}.c,3));' ...
                    'clear tmpdat tmpind;' ];
    cb_changeradii = [  'tmpdat = get(gcbf, ''userdata'');' ...
                       'tmpdat.vol.r = str2num(get(gcbo, ''string''));' ...
                       'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeconduct = [  'tmpdat = get(gcbf, ''userdata'');' ...
                       'tmpdat.vol.c = str2num(get(gcbo, ''string''));' ...
                       'set(gcf, ''userdata'', tmpdat)' ];
    cb_changeorigin = [  'tmpdat = get(gcbf, ''userdata'');' ...
                       'tmpdat.vol.o = str2num(get(gcbo, ''string''));' ...
                       'set(gcf, ''userdata'', tmpdat)' ];
    % cb_fitelec = [ 'if get(gcbo, ''value''),' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''off'');' ...
    %                'else' ...
    %                '  set(findobj(gcbf, ''tag'', ''origin''), ''enable'', ''on'');' ...
    %                'end;' ];
    
    userdata    = [];
    
    geomvert = [1 1 1 1 1 1 1];
    
    geomhorz = {
        [1]
        [1 1.3]
        [1 1.3]
        [1.28 1 0.6]
        [1]
        [1]
        [1]
               };
    
    % define each individual graphical user element
    elements  = { ...
        { 'style' 'text'        'string' 'Volume conductor model (brain CVF skull skin):' }  ...
        { 'style' 'text'        'string' 'Shell radii (mm)' }  ...
        { 'style' 'edit'        'string'  num2str(defaultvolume.r,3) 'tag' 'radii' 'callback', cb_changeradii} ...
        { 'style' 'text'        'string' 'Shell conductances (S/m)' }  ...
        { 'style' 'edit'        'string'  num2str(defaultvolume.c,3) 'tag' 'conduct' 'callback', cb_changeconduct} ...
        { 'style' 'text'        'string' 'Omit channels for dipole fit' } ...
        { 'style' 'edit'        'string' ''  'tag' 'elec' } ...
        { 'style' 'pushbutton'  'string' 'List' 'callback' cb_selectelectrodes } ... 
        { } ...
        { 'style' 'text'        'string' 'Note: Under menu item ''Edit > Channel locations,'' use the ''3-D Center''' }  ...
        { 'style' 'text'        'string' 'button to fit the center of the spherical head model before fitting dipoles.' } ...
                };
    
    result = inputgui( geomhorz, elements, 'pophelp(''pop_dipfit_settings'')', ...
                                     'Dipole fit settings - pop_dipfit_settings()', userdata, 'normal', geomvert );
    
    if isempty(result), return; end
    options = {};
    options = { options{:} 'radii'        str2num(result{1}) };
    options = { options{:} 'conductances' str2num(result{2}) };
    options = { options{:} 'electrodes'   setdiff(1:EEG.nbchan, str2num(result{3})) };

else
    options = varargin;
end

g = finputcheck( options, { 'radii'        'float'     []             [defaultvolume.r];
                            'conductances' 'float'     []             [defaultvolume.c];
                            'origin'       'float'     []             [defaultvolume.o];
                            'projectelectrodes' 'integer' [0 1]       [0];
                            'electrodes'   'integer'   [1 Inf]        [1:EEG.nbchan]});
if isstr(g), error(g); end;

% remember the volume conductor and electrode settings
OUTEEG.dipfit.chansel = g.electrodes;
OUTEEG.dipfit.vol.r   = g.radii;
OUTEEG.dipfit.vol.c   = g.conductances;
OUTEEG.dipfit.vol.o   = defaultvolume.o;

if g.projectelectrodes & isfield(EEG, 'chanlocs')
  if isfield(OUTEEG.chanlocs, 'X') & isfield(OUTEEG.chanlocs, 'Y') & isfield(OUTEEG.chanlocs, 'Z')
    % here we have the option to remember the original channel locations
    % but it is not decided yet upon where to put them
    % OUTEEG.urchanlocs = OUTEEG.chanlocs;

    % project the electrodes towards the skin
    for el = 1:length(EEG.chanlocs)
      xyz(1) = OUTEEG.chanlocs(el).X;
      xyz(2) = OUTEEG.chanlocs(el).Y;
      xyz(3) = OUTEEG.chanlocs(el).Z;
      xyz = xyz - g.origin;
      xyz = max(g.radii) * xyz / norm(xyz);
      xyz = xyz + g.origin;
      OUTEEG.chanlocs(el).X = xyz(1);
      OUTEEG.chanlocs(el).Y = xyz(2);
      OUTEEG.chanlocs(el).Z = xyz(3);
    end
    OUTEEG.chanlocs = convertlocs(OUTEEG.chanlocs, 'cart2all');
  else
    error('could not project electrodes because carthesian coordinates are not available');
  end
else 
    warning off backtrace; % avoid long error messages
end

com = sprintf('%s = pop_dipfit_settings( %s, %s);', inputname(1), inputname(1), vararg2str(options));
