% egg_multifit - fit multiple component dipoles using DIPFIT 
%
% Usage:
%         >> EEG = eeg_multifit(EEG); % fit all dipoles
%         >> EEG = eeg_multifit(EEG, 'key', 'val', ...);
%
% Inputs:
%  EEG      - input dataset.
%
% Optional inputs:
%  'comps'    - [integer array] fit selected components only.
%  'dipfit'   - [cell array] dipfit settings. Default is 
%               { 'projectelectrodes', 1 }. Note that this option
%               currently might modify the coordinates of the channels.
%
% Outputs:
%  EEG      - output dataset with updated "EEG.dipfit" field
%
% Note: residual variance is set to NaN if Dipfit crashes
%  
% Author: Arnaud Delorme, SCCN/INC/UCSD, La Jolla, Sept. 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 9/2003 Arnaud Delorme, SCCN/INC/UCSD, arno@sccn.ucsd.edu
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

% $Log: eeg_multifit.m,v $
% Revision 1.8  2005/05/24 18:00:14  arno
% remove cell2mat
%
% Revision 1.7  2003/11/07 00:30:01  arno
% remving grid
%
% Revision 1.6  2003/11/05 16:22:03  arno
% homogenous -> homogeneous
%
% Revision 1.5  2003/10/07 01:14:14  arno
% *** empty log message ***
%
% Revision 1.4  2003/10/04 02:02:24  arno
% code for removing the eye channels
%
% Revision 1.3  2003/10/04 01:46:30  arno
% error if trying to fit more than 69 electrodes
%
% Revision 1.2  2003/10/03 23:36:15  arno
% typo
%
% Revision 1.1  2003/09/29 23:06:08  arno
% Initial revision
%

function EEG = eeg_multifit(EEG, varargin);
    
    if nargin < 1
        help eeg_multifit;
        return;
    end;
    
    % checking parameters
    % -------------------
    ncomps = size(EEG.icaweights,1);
    if ncomps == 0, error('you must run ICA first'); end;
    g = finputcheck(varargin, { 'comps'   'integer'  [1 ncomps]   [1:ncomps];
                                'dipfit'  'cell'     []        {'projectelectrodes', 1} });
    
    if isstr(g), error(g); end;
    
    %if isfield(EEG, 'dipfit')
    %    EEG     = rmfield(EEG, 'dipfit');
    %end;
    EEG     = eeg_checkset(EEG, 'chanlocs_homogeneous');
    if length(EEG.chanlocs) > 69 & length(g.dipfit) < 3
        [tmp tmp2 indreye] = intersect( 'reye', lower({ EEG.chanlocs.labels }));
        [tmp tmp2 indleye] = intersect( 'leye', lower({ EEG.chanlocs.labels }));
        elecsel = setdiff([1:71], [indleye indreye]);
        disp('Removing eye channels');
        EEG     = pop_dipfit_settings( EEG, g.dipfit{:}, 'electrodes', elecsel);
    else
        EEG     = pop_dipfit_settings( EEG, g.dipfit{:});        
    end;
    if length(EEG.dipfit.chansel) > 69, error('Cannot fit more than 69 electrodes'); end;
    
    EEG     = pop_dipfit_batch( EEG, g.comps); %, linspace(-1,1,10), linspace(-1,1,10), linspace(-0.5,1,10));
    %fprintf('Local fitting:');
    %for index = g.comps
    %    fprintf('%d\t', index);
    %    aroundx = EEG.dipfit.model(index).posxyz(1,1)+linspace(-0.1, 0.1, 10);
    %    aroundy = EEG.dipfit.model(index).posxyz(1,2)+linspace(-0.1, 0.1, 10);
    %    aroundz = EEG.dipfit.model(index).posxyz(1,3)+linspace(-0.1, 0.1, 10);
    %    EEG     = pop_dipfit_batch( EEG, index, aroundx, aroundy, aroundz);
    %end;
    %fprintf('\n');    
    %return;
    
    chansel =  EEG.dipfit.chansel;
    elc     = [ [EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]' ];
    for i = g.comps(:)'
        EEG.dipfit.model(i).active = [1];
        EEG.dipfit.model(i).select = [1];
        try,
            EEG.dipfit.model(i) = dipfit_manual(EEG.dipfit.model(i), EEG.icawinv(chansel,i), ...
                                                elc(chansel,:), EEG.dipfit.vol, 'method', 'position' );
        catch, EEG.dipfit.model(i).rv =NaN; 
        end;
    end;
    return;