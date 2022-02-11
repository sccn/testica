% geticadata - get data from subject and associated ICA decomposition
%
% Usage:
%   EEG = geticadata( dataset_num, algostr );
%
% Inputs:
%   dataset_num  - [integer] dataset number from 1 to 14
%   algostr      - [string] ICA algorithm name
%
% Output:
%   EEG          - EEGLAB dataset
%
% Notes:
% possible algorithms are 
% algorithms = { 'runica ext' 'runica' 'erica' 'sons' 'shibbs' 'jade' ...
%                'tica' 'nat. grad.' 'jade opt.' 'jade td' 'eea' 'eobi' ...
%                'acsobiro' 'evd 24' 'evd' 'sobi' 'amuse' 'pearson' ...
%                'icaml' 'fastica' 'pca' 'icams' 'simbec' };
%
% Author: Arnaud Delorme, SCCN, 2007

function EEG = geticadata(dat, algoname, TMPEEG);
    
% DATASET IS SET so dataset is read below
if nargin < 3
    DATASET = dat;
    processdat;
else
    processdat;
    EEG = TMPEEG;
end;

% load data for specific ICA algorithm

tmpind = cellfun('isempty', { allalgs.algo });
for ind = find(tmpind), allalgs(ind).name  = ''; allalgs(ind).name  = ''; end;
ALGONUM = strmatch(algoname, { allalgs.name }, 'exact');
if isempty(ALGONUM), ALGONUM = strmatch(algoname, { allalgs.algo }, 'exact'); end;
if length(ALGONUM) ~= 1, error('Algorithm not found'); end;
disp('*****************************************************');
filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',dat, ALGONUM, allalgs(ALGONUM).algo);
disp(['Loading ICA weight file ' filename ]);
disp('*****************************************************');
load('-mat', filename);
EEG.icaweights = W;
if size(W,1) ~= size(W,2), error('Non square weight matrix'); end;
EEG.icasphere  = eye(size(W));
EEG.icawinv    = pinv(W);
EEG.icaact     = [];
EEG.subject    = algoname;
EEG.setname    = [EEG.setname ' ' algoname ];
for i = 1:length(allrv)
    EEG.dipfit.model(i).rv = allrv(i);
    EEG.dipfit.model(i).posxyz = allposxyz(:,i)';
    EEG.dipfit.model(i).momxyz = allmomxyz(:,i)';
end;
EEG = eeg_checkset(EEG);




