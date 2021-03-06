% pop_dipfit_manual() - interactively do dipole fit of selected ICA components
%
% Usage: 
%  >> OUTEEG = pop_dipfit_manual( INEEG )
%
% Inputs:
%   INEEG       input dataset
%
% Outputs:
%   OUTEEG      output dataset
%
% Author: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%         Arnaud Delorme, SCCN, La Jolla 2003

% SMI, University Aalborg, Denmark http://www.smi.auc.dk/
% FC Donders Centre, University Nijmegen, the Netherlands http://www.fcdonders.kun.nl/

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

% $Log: pop_dipfit_manual.m,v $
% Revision 1.29  2005/05/24 18:03:11  arno
% remove cell2mat
%
% Revision 1.28  2004/08/12 22:03:40  arno
% same
%
% Revision 1.27  2004/08/12 22:02:20  arno
% same
%
% Revision 1.26  2004/08/12 22:00:45  arno
% topoplot with dipole
%
% Revision 1.25  2004/08/12 21:52:59  arno
% do not crash
%
% Revision 1.24  2004/03/26 01:33:29  arno
% plot only active dipoles
%
% Revision 1.23  2003/12/04 17:09:07  arno
% num2str -> str2num
%
% Revision 1.22  2003/12/04 16:05:58  arno
% setting values when modifying gui
%
% Revision 1.21  2003/10/31 16:45:52  arno
% worind from Scott's feedback
%
% Revision 1.20  2003/10/30 02:39:14  arno
% gui typo
%
% Revision 1.19  2003/10/29 23:21:01  arno
% more checking for component index
%
% Revision 1.18  2003/10/29 23:03:37  arno
% [Adecrease window width
%
% Revision 1.17  2003/10/29 16:08:00  arno
% removing debug msg
%
% Revision 1.16  2003/10/29 03:12:11  arno
% contrain electrode to sphere
%
% Revision 1.15  2003/09/30 15:42:34  roberto
% minor bug fix, related to the change in dipole constraint handling
%
% Revision 1.14  2003/09/12 08:43:18  roberto
% changed symmetry constraint into optional input argument
%
% Revision 1.12  2003/07/01 23:49:08  arno
% implementing dipole flipping
%
% Revision 1.11  2003/06/30 02:11:54  arno
% *** empty log message ***
%
% Revision 1.10  2003/06/16 10:10:11  roberto
% added interruptible dialog (gui) for non-linear fit
%
% Revision 1.9  2003/06/13 16:48:57  arno
% remove chanlocs conversion
%
% Revision 1.8  2003/06/13 01:26:01  arno
% automatic conversion of channel location files
%
% Revision 1.7  2003/06/13 01:17:49  arno
% adding plotting button
%
% Revision 1.6  2003/03/12 10:32:32  roberto
% fixed dialog title
%
% Revision 1.5  2003/03/06 15:58:21  roberto
% fixed bug with channel selection of EEG data
%
% Revision 1.3  2003/03/03 16:52:27  roberto
% modified for posxyz/momxyz instead of dip.pos/dip.mom
% changed large listbox into edit field
%
% Revision 1.2  2003/02/28 23:03:14  arno
% making the listbox to scroll component
%
% Revision 1.1  2003/02/24 10:06:15  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_manual( EEG, subfunction, parent, dipnum )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% the code for this interactive dialog has 4 major parts
% - draw the graphical user interface
% - synchronize the gui with the data
% - synchronize the data with the gui
% - execute the actual dipole analysis

if nargin<1
  help pop_dipfit_manual;
  return
elseif nargin==1

OUTEEG = EEG;
com = '';

if ~isfield(EEG, 'chanlocs')
  error('No electrodes present');
end

if ~isfield(EEG, 'icawinv')
  error('No ICA components to fit');
end

if ~isfield(EEG, 'dipfit')
  error('General dipolefit settings not specified');
end

if ~isfield(EEG.dipfit, 'vol')
  error('Dipolefit volume conductor model not specified');
end

% select all ICA components as 'fitable'
select = 1:size(EEG.icawinv,2);
if ~isfield(EEG.dipfit, 'current')
  EEG.dipfit.current = 1;
end

% verify the presence of a dipole model
if ~isfield(EEG.dipfit, 'model')
  % create empty dipole model for each component
  for i=select
    EEG.dipfit.model(i).posxyz = zeros(2,3);
    EEG.dipfit.model(i).momxyz = zeros(2,3);
    EEG.dipfit.model(i).rv     = 1;
    EEG.dipfit.model(i).active = [1];
    EEG.dipfit.model(i).select = [1];
  end
end

% verify the size of each dipole model
for i=select
  if ~isfield(EEG.dipfit.model, 'posxyz') | length(EEG.dipfit.model) < i | isempty(EEG.dipfit.model(i).posxyz)
    % replace all empty dipole models with a two dipole model, of which one is active
    EEG.dipfit.model(i).active = [1];
    EEG.dipfit.model(i).select = [1];
    EEG.dipfit.model(i).rv = 1;
    EEG.dipfit.model(i).posxyz = zeros(2,3);
    EEG.dipfit.model(i).momxyz = zeros(2,3);
  elseif size(EEG.dipfit.model(i).posxyz,1)==1
    % replace all one dipole models with a two dipole model
    EEG.dipfit.model(i).active = [1];
    EEG.dipfit.model(i).select = [1];
    EEG.dipfit.model(i).posxyz = [EEG.dipfit.model(i).posxyz; [0 0 0]];
    EEG.dipfit.model(i).momxyz = [EEG.dipfit.model(i).momxyz; [0 0 0]];
  elseif size(EEG.dipfit.model(i).posxyz,1)>2
    % replace all more-than-two dipole models with a two dipole model
    warning('pruning dipole model to two dipoles');
    EEG.dipfit.model(i).active = [1];
    EEG.dipfit.model(i).select = [1];
    EEG.dipfit.model(i).posxyz = EEG.dipfit.model(i).posxyz(1:2,:);
    EEG.dipfit.model(i).momxyz = EEG.dipfit.model(i).momxyz(1:2,:);
  end
end

% default is not to use symmetry constraint
constr = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct the graphical user interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the callback functions for the interface elements
cb_plotmap         = 'pop_dipfit_manual(EEG, ''dialog_plotmap'', gcbf);';
cb_selectcomponent = 'pop_dipfit_manual(EEG, ''dialog_selectcomponent'', gcbf);';
cb_checkinput      = 'pop_dipfit_manual(EEG, ''dialog_checkinput'', gcbf);';
cb_fitposition     = 'pop_dipfit_manual(EEG, ''dialog_getvalue'', gcbf); pop_dipfit_manual(EEG, ''dipfit_position'', gcbf); pop_dipfit_manual(EEG, ''dialog_setvalue'', gcbf);';
cb_fitmoment       = 'pop_dipfit_manual(EEG, ''dialog_getvalue'', gcbf); pop_dipfit_manual(EEG, ''dipfit_moment''  , gcbf); pop_dipfit_manual(EEG, ''dialog_setvalue'', gcbf);';
cb_close           = 'close(gcbf)';
cb_help            = 'pophelp(''pop_dipfit_manual'');';
cb_ok              = 'uiresume(gcbf);'; 
cb_plotdip         = 'pop_dipfit_manual(EEG, ''dialog_getvalue'', gcbf); pop_dipfit_manual(EEG, ''dialog_plotcomponent'', gcbf);';
cb_flip1           = 'pop_dipfit_manual(EEG, ''dialog_flip'', gcbf, 1);';
cb_flip2           = 'pop_dipfit_manual(EEG, ''dialog_flip'', gcbf, 2);';

% vertical layout for each line 
geomvert  =   [1 1 1 1 1 1 1 1 1]; 

% horizontal layout for each line 
geomhoriz = {
	      [0.8 0.5 0.8 1 1]
              [1]
              [0.7 0.7  0.7     2 2 1]
              [0.7 0.5 0.2 0.5 0.2 2 2 1]
              [0.7 0.5 0.2 0.5 0.2 2 2 1]
              [1] 
              [1 1 1]
              [1]
              [1 1 1]
            };

% define each individual graphical user element
elements  = { ...
              { 'style' 'text'       'string' 'Component to fit'    } ...
              { 'style' 'edit'       'string' 'dummy'    'tag' 'component' 'callback' cb_selectcomponent } ...
              { 'style' 'pushbutton' 'string' 'Plot map'                   'callback' cb_plotmap } ...
              { 'style' 'text'       'string' 'Residual variance = '                                             } ...
              { 'style' 'text'       'string' 'dummy'         'tag' 'relvar'                                  } ...
              { } ...
              { 'style' 'text'    'string' 'dipole'     } ...
              { 'style' 'text'    'string' 'active'     } ...
              { 'style' 'text'    'string' 'fit'        } ...
              { 'style' 'text'    'string' 'position'   } ...
              { 'style' 'text'    'string' 'moment'     } ...
              { } ...
              ...
              { 'style' 'text'        'string' '#1' 'tag' 'dip1'                                 } ...
              { 'style' 'checkbox'    'string' ''   'tag' 'dip1act'    'callback' cb_checkinput  } { } ...
              { 'style' 'checkbox'    'string' ''   'tag' 'dip1sel'    'callback' cb_checkinput  } { } ...
              { 'style' 'edit'        'string' ''   'tag' 'dip1pos'    'callback' cb_checkinput  } ...
              { 'style' 'edit'        'string' ''   'tag' 'dip1mom'    'callback' cb_checkinput  } ...
              { 'style' 'pushbutton'  'string' 'Flip (in|out)'         'callback' cb_flip1       } ...   
              ...
              { 'style' 'text'        'string' '#2' 'tag' 'dip2'                                 } ...
              { 'style' 'checkbox'    'string' ''   'tag' 'dip2act'    'callback' cb_checkinput  } { } ...
              { 'style' 'checkbox'    'string' ''   'tag' 'dip2sel'    'callback' cb_checkinput  } { } ...
              { 'style' 'edit'        'string' ''   'tag' 'dip2pos'    'callback' cb_checkinput  } ...
              { 'style' 'edit'        'string' ''   'tag' 'dip2mom'    'callback' cb_checkinput  } ...
              { 'style' 'pushbutton'  'string' 'Flip (in|out)'         'callback' cb_flip2       } ...   
              ...
              { } { 'style' 'checkbox' 'string' 'Symmetry constrain for dipole #2' 'tag' 'dip2sym' 'callback' cb_checkinput  'value' 1 } ...
              { } { } { } ...
              { 'style' 'pushbutton'  'string' 'Fit selected dipole position(s)' 'callback' cb_fitposition } ...   
              { 'style' 'pushbutton'  'string' 'Fit selected dipole moment(s)'   'callback' cb_fitmoment   } ...   
              { 'style' 'pushbutton'  'string' 'Plot dipole(s)'                  'callback' cb_plotdip     } ...   
	    };

% add the cancel, help and ok buttons at the bottom

geomvert  = [geomvert 1 1];

geomhoriz = {geomhoriz{:} [1] [1 1 1]};

elements  = { elements{:} ...
              { } ...
              { 'Style', 'pushbutton', 'string', 'Cancel', 'callback', cb_close } ...
              { 'Style', 'pushbutton', 'string', 'Help',   'callback', cb_help  } ...
              { 'Style', 'pushbutton', 'string', 'OK',     'callback', cb_ok    } ...
            };

% activate the graphical interface
supergui(0, geomhoriz, geomvert, elements{:});
dlg = gcf;
set(gcf, 'name', 'Manual dipole fit -- pop_dipfit_manual()');
set(gcf, 'userdata', EEG);
pop_dipfit_manual(EEG, 'dialog_setvalue', dlg);
uiwait(dlg);
if ishandle(dlg)
  pop_dipfit_manual(EEG, 'dialog_getvalue', dlg);
  % FIXME, rv is undefined since the user may have changed dipole parameters
  % FIXME, see also dialog_getvalue subfucntion
  OUTEEG = get(dlg, 'userdata');
  close(dlg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implement all subfunctions through a switch-yard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif nargin>=3
  %disp(subfunction)
  EEG = get(parent, 'userdata');

  switch subfunction

  case 'dialog_selectcomponent'
    current = get(findobj(parent, 'tag', 'component'), 'string');
    current = str2num(current);
    current = current(1);
    current = min(current, size(EEG.icaweights,1));
    current = max(current, 1);
    set(findobj(parent, 'tag', 'component'), 'string', int2str(current));
    EEG.dipfit.current = current;
    % reassign the global EEG object back to the dialogs userdata
    set(parent, 'userdata', EEG);
    % redraw the dialog with the current model
    pop_dipfit_manual(EEG, 'dialog_setvalue', parent);
    
   case 'dialog_plotmap'
    current = str2num(get(findobj(parent, 'tag', 'component'), 'string'));
    figure; pop_topoplot(EEG, 0, current, [ 'IC ' num2str(current) ], [1 1], 1); 

   case 'dialog_plotcomponent'
    current = get(findobj(parent, 'tag', 'component'), 'string');
    EEG.dipfit.current = str2num(current);
    if ~isempty( EEG.dipfit.current )
        pop_dipplot(EEG, 'DIPFIT',  EEG.dipfit.current, 'normlen', 'on', 'image', 'mri');
    end;
    
   case 'dialog_checkinput'
    if get(findobj(parent, 'tag', 'dip1sel'), 'value') & ~get(findobj(parent, 'tag', 'dip1act'), 'value')
        set(findobj(parent, 'tag', 'dip1act'), 'value', 1);
    end
    if get(findobj(parent, 'tag', 'dip2sel'), 'value') & ~get(findobj(parent, 'tag', 'dip2act'), 'value')
        set(findobj(parent, 'tag', 'dip2act'), 'value', 1);
    end
    if ~all(size(str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip1pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:)));
    else
        EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string'));
    end
    if ~all(size(str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip2pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:)));
    else
        EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string'));
    end
    if ~all(size(str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip1mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:)));
    else
        EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string'));
    end
    if ~all(size(str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string')))==[1 3])
        set(findobj(parent, 'tag', 'dip2mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:)));
    else
        EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string'));
    end
    if get(findobj(parent, 'tag', 'dip2sel'), 'value') & get(findobj(parent, 'tag', 'dip2sym'), 'value') & ~get(findobj(parent, 'tag', 'dip1sel'), 'value')
        set(findobj(parent, 'tag', 'dip2sel'), 'value', 0);
    end
    set(parent, 'userdata', EEG);

  case 'dialog_setvalue'
    % synchronize the gui with the data
    set(findobj(parent, 'tag', 'component'), 'string', EEG.dipfit.current);
    set(findobj(parent, 'tag', 'relvar' ), 'string', sprintf('%0.2f%%', EEG.dipfit.model(EEG.dipfit.current).rv * 100));
    set(findobj(parent, 'tag', 'dip1act'), 'value', ismember(1, EEG.dipfit.model(EEG.dipfit.current).active));
    set(findobj(parent, 'tag', 'dip2act'), 'value', ismember(2, EEG.dipfit.model(EEG.dipfit.current).active));
    set(findobj(parent, 'tag', 'dip1sel'), 'value', ismember(1, EEG.dipfit.model(EEG.dipfit.current).select));
    set(findobj(parent, 'tag', 'dip2sel'), 'value', ismember(2, EEG.dipfit.model(EEG.dipfit.current).select));
    set(findobj(parent, 'tag', 'dip1pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(1,:)));
    set(findobj(parent, 'tag', 'dip2pos'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).posxyz(2,:)));
    set(findobj(parent, 'tag', 'dip1mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(1,:)));
    set(findobj(parent, 'tag', 'dip2mom'), 'string', sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(2,:)));

  case 'dialog_getvalue'
    % synchronize the data with the gui
    if get(findobj(parent, 'tag', 'dip1act'), 'value'); active = [1]; else active = []; end; 
    if get(findobj(parent, 'tag', 'dip2act'), 'value'); active = [active 2]; end; 
    if get(findobj(parent, 'tag', 'dip1sel'), 'value'); select = [1]; else select = []; end; 
    if get(findobj(parent, 'tag', 'dip2sel'), 'value'); select = [select 2]; end; 
    posxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1pos'), 'string')); 
    posxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2pos'), 'string')); 
    momxyz(1,:) = str2num(get(findobj(parent, 'tag', 'dip1mom'), 'string')); 
    momxyz(2,:) = str2num(get(findobj(parent, 'tag', 'dip2mom'), 'string')); 
    % assign the local values to the global EEG object
    EEG.dipfit.model(EEG.dipfit.current).posxyz = posxyz;
    EEG.dipfit.model(EEG.dipfit.current).momxyz = momxyz;
    EEG.dipfit.model(EEG.dipfit.current).active = active;
    EEG.dipfit.model(EEG.dipfit.current).select = select;
    % FIXME, rv is undefined after a manual change of parameters
    % FIXME, this should either be undated continuously or upon OK buttonpress
    % EEG.dipfit.model(EEG.dipfit.current).rv = nan;

    % reassign the global EEG object back to the dialogs userdata
    set(parent, 'userdata', EEG);

  case 'dialog_flip'
    % make a copy of the electrodes, with only the selected ones
    current = EEG.dipfit.current;
    moment  = EEG.dipfit.model(current).momxyz;
    EEG.dipfit.model(current).momxyz(dipnum,:) = [ -moment(dipnum,1) -moment(dipnum,2) -moment(dipnum,3)];
    set(findobj(parent, 'tag', ['dip' int2str(dipnum) 'mom']), 'string', ...
                      sprintf('%0.3f %0.3f %0.3f', EEG.dipfit.model(EEG.dipfit.current).momxyz(dipnum,:)));
    set(parent, 'userdata', EEG);
   
  case 'dipfit_position'
    % make a copy of the electrodes and get the selected channels
    elc     = getelecpos(EEG.chanlocs, EEG.dipfit);
    chansel = EEG.dipfit.chansel;
    elc     = elc(chansel,:);
    % make some other local copies for convenience
    current = EEG.dipfit.current;
    model   = EEG.dipfit.model(current);
    pot     = EEG.icawinv(chansel,current);
    vol     = EEG.dipfit.vol;
    if isempty(model.select)
      % no dipoles are selected, nothing to do
      return;
    end
    % set the symmetry constraint for dipole 2
    dipfitdefs;
    if get(findobj(parent, 'tag', 'dip2sym'), 'value') & get(findobj(parent, 'tag', 'dip2sel'), 'value')
      constr = defaultconstraint;
    else
      constr = [];
    end
    % make a dialog which can interrupt the fitting procedure
    fig = figure('visible', 'off');
    supergui( fig, {1 1}, [], ...
        {'style' 'text' 'string' 'Press button below to stop fitting' }, ...
        {'style' 'pushbutton' 'string' 'Interupt' 'callback' 'figure(gcbf); set(gcbf, ''tag'', ''stop'');' } );
    drawnow;
    try, 
        model = dipfit_manual(model, pot, elc, vol, 'method', 'position', 'constraint', constr );
        % remember the fitted dipole object
        EEG.dipfit.model(current) = model;
        % reassign the global EEG object back to the dialogs userdata
        set(parent, 'userdata', EEG);
    catch,
        disp('Dipole localization failed');
    end;
    if ishandle(fig)
      close(fig);
    end

  case 'dipfit_moment'
    % make a copy of the electrodes, with only the selected ones
    elc     = getelecpos(EEG.chanlocs, EEG.dipfit);
    chansel = EEG.dipfit.chansel;
    elc     = elc(chansel,:);
    % make some other local copies for convenience
    current = EEG.dipfit.current;
    model   = EEG.dipfit.model(current);
    pot     = EEG.icawinv(chansel,current);
    vol     = EEG.dipfit.vol;
    if isempty(model.select)
      % no dipoles are selected, nothing to do
      return;
    end
    model = dipfit_manual(model, pot, elc, vol, 'method', 'moment' );
    % remember the fitted dipole object
    EEG.dipfit.model(current) = model;
    % reassign the global EEG object back to the dialogs userdata
    set(parent, 'userdata', EEG);

  otherwise
    error('unknown subfunction for pop_dipfit_manual');
  end		% switch subfunction

end		% if nargin

% get electrode positions from eeglag
% -----------------------------------
function elc = getelecpos(chanlocs, dipfitstruct);
    try,
        elc =  [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]' ];
    catch
        disp('No 3-D carthesian coordinates; re-computing them from 2-D polar coordinates');
        EEG.chanlocs = convertlocs(EEG.chanlocs, 'topo2all');
        elc =  [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]' ];
    end;
    % constrain electrode to sphere
    % -----------------------------
    disp('Constraining electrodes to sphere');
    elc = elc - repmat( dipfitstruct.vol.o, [size(elc,1) 1]); % recenter 
    % (note the step above is not needed since the origin should always be 0)
    elc = elc ./ repmat( sqrt(sum(elc.*elc,2)), [1 3]); % normalize
    elc = elc * max(dipfitstruct.vol.r);         % head size
    
    %for index= 1:size(elc,1)
    %    elc(index,:) = max(dipfitstruct.vol.r) * elc(index,:) /norm(elc(index,:));
    %end;
