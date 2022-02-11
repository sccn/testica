% pop_dipfit_batch() - interactively do batch scan of all ICA components
%                      with a single dipole
%
% Usage: 
%  >> OUTEEG = pop_dipfit_batch( INEEG ); % pop up interactive window
%  >> OUTEEG = pop_dipfit_batch( INEEG, comps );
%  >> OUTEEG = pop_dipfit_batch( INEEG, comps, xgrid, ygrid, zgrid, thresh )
%
% Inputs:
%   INEEG     - input dataset
%   comps     - [integer array] component indices
%   xgrid     - [float array] x-grid. Default is 10 elements between
%               -1 and 1.
%   ygrid     - [float array] y-grid. Default is 10 elements between
%               -1 and 1.
%   zgrid     - [float array] z-grid. Default is 10 elements between
%               -1 and 1.
%   threshold - [float] threshold in percent. Default 40.
%
% Outputs:
%   OUTEEG      output dataset
%
% Authors: Robert Oostenveld, SMI/FCDC, Nijmegen 2003
%          Arnaud Delorme, SCCN, La Jolla 2003

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

% $Log: pop_dipfit_batch.m,v $
% Revision 1.17  2005/05/24 17:59:24  arno
% remove cell2mat
%
% Revision 1.16  2005/03/17 00:49:42  arno
% debug for electrode with no locations
%
% Revision 1.15  2003/11/25 01:06:30  arno
% header edit
%
% Revision 1.14  2003/11/25 00:58:58  arno
% interface text
%
% Revision 1.13  2003/10/29 22:44:44  arno
% wording of GUI
%
% Revision 1.12  2003/10/29 16:39:27  arno
% fixing return command
%
% Revision 1.11  2003/10/29 03:13:32  arno
% constrain electrodes to sphere
%
% Revision 1.10  2003/10/08 17:35:19  arno
% fixing thresh command line bug
%
% Revision 1.9  2003/10/07 00:45:03  arno
% typo
%
% Revision 1.8  2003/09/29 22:17:07  arno
% debuging new select option
%
% Revision 1.7  2003/08/05 17:50:14  arno
% fixing dipole index for batch
%
% Revision 1.6  2003/07/01 22:23:01  arno
% nothing
%
% Revision 1.5  2003/06/30 01:46:30  arno
% fixing waitbar
%
% Revision 1.4  2003/06/30 01:43:54  arno
% *** empty log message ***
%
% Revision 1.3  2003/03/06 15:58:02  roberto
% fixed bug with channel selection of EEG data
% changed rejection threshold into percent
%
% Revision 1.1  2003/02/24 10:05:58  roberto
% Initial revision
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OUTEEG, com] = pop_dipfit_batch(EEG, select, xgrid, ygrid, zgrid, reject );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help pop_dipfit_batch;
   return;
end;

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

dipfitdefs
if nargin < 2
    % get the default values and filenames
    promptstr = { 'Component(s) (not faster if few comp.)', ...
                  'Grid in X-direction', ...
                  'Grid in Y-direction', ...
                  'Grid in Z-direction', ...
                  'Rejection threshold RV(%)' };
    
    inistr    = { 
        [ '1:' int2str(size(EEG.icawinv,2)) ], ...
        xgridstr, ...
        ygridstr, ...
        zgridstr, ...
        rejectstr };

    result = inputdlg2( promptstr, 'Batch dipole fit -- pop_dipfit_batch()', 1,  inistr, 'pop_dipfit_batch');
    
    if length(result)==0
        % user pressed cancel
        return
    end
    
    select = eval( [ '[' result{1} ']' ]);
    xgrid  = eval( result{2} );
    ygrid  = eval( result{3} );
    zgrid  = eval( result{4} );
    reject = eval( result{5} ) / 100;	% string is in percent
    options = { };
else
    if nargin < 2
        select = [1:size(EEG.icawinv,2)];
    end;
    if nargin < 3
        xgrid  = eval( xgridstr );
    end;
    if nargin < 4
        ygrid  = eval( ygridstr );
    end;
    if nargin < 5
        zgrid  = eval( zgridstr );
    end;
    if nargin < 6
        reject  = eval( rejectstr );
    end;
    options = { 'waitbar' 'none' };
end;
    
% make a copy of all electrodes
% -----------------------------
elc = getelecpos(EEG.chanlocs, EEG.dipfit);
chansel = EEG.dipfit.chansel;

% perform batch fit with single dipole for all selected channels and components
model  = dipfit_batch(EEG.icawinv(chansel,select), elc(chansel,:), EEG.dipfit.vol, ...
                      'xgrid', xgrid, 'ygrid', ygrid, 'zgrid', zgrid, options{:});

% discard the models that do not have sufficient quality with a single dipole
model  = dipfit_reject(model, reject);

% subsequent scanning with a symmetric dipole pair could automatically 
% be done on those components that have a poor fit with a single dipole

for index = 1:length(model)
    model(index).active = 1;
    model(index).select = 1;
end;

OUTEEG.dipfit.model(select) = model;

com = sprintf('%s = pop_dipfit_batch( %s,%s);', inputname(1), inputname(1), vararg2str({ select xgrid ygrid zgrid reject }));

% get electrode positions from eeglag
% -----------------------------------
function elc = getelecpos(chanlocs, dipfitstruct);
    for index = 1:length(chanlocs)
        if isempty(chanlocs(index).X)
            chanlocs(index).X = NaN;
            chanlocs(index).Y = NaN;
            chanlocs(index).Z = NaN;
        end;
    end;
            
    elc     = [ [chanlocs.X]' [chanlocs.Y]' [chanlocs.Z]' ];

    % constrain electrode to sphere
    % -----------------------------
    disp('Constraining electrodes to sphere');
    elc = elc - repmat( dipfitstruct.vol.o, [size(elc,1) 1]); % recenter 
    % (note the step above is not needed since the origin should always be 0)
    elc = elc ./ repmat( sqrt(sum(elc.*elc,2)), [1 3]); % normalize
    elc = elc * max(dipfitstruct.vol.r);         % head size
    