% have to define 2 values
% DATASET: sample dataset to use
% ALGONUM: algorythm number

if ~exist('pop_loadset.m')
    eeglab; close;
end;
clear ALLEEG EEG CURRENTSET ALLCOM LASTCOM;

parentpath = fileparts(which('processdat'));
parentpath = fullfile(parentpath, 'datasets');
if exist('DATASET') == 1
    switch DATASET
     case 1 , EEG = pop_loadset( 'km81.set', parentpath );
     case 2 , EEG = pop_loadset( 'jo74.set', parentpath );
     case 3 , EEG = pop_loadset( 'ds76.set', parentpath );
     case 4 , EEG = pop_loadset( 'cj82.set', parentpath );
     case 5 , EEG = pop_loadset( 'ap82.set', parentpath );
     case 6 , EEG = pop_loadset( 'ke70.set', parentpath );
     case 7 , EEG = pop_loadset( 'tp62.set', parentpath );
     case 8 , EEG = pop_loadset( 'cz84.set', parentpath );
     case 9 , EEG = pop_loadset( 'gm84.set', parentpath ); % good to test extended
     case 10, EEG = pop_loadset( 'gv84.set', parentpath );
     case 11, EEG = pop_loadset( 'nf68.set', parentpath );
     case 12, EEG = pop_loadset( 'ds80.set', parentpath );
     case 13, EEG = pop_loadset( 'kb77.set', parentpath );
     case 14, EEG = pop_loadset( 'ts79.set', parentpath );
    end;
    EEG.icaweights = [];
    EEG.icasphere  = [];
    EEG.icawinv    = [];
    EEG.icaact     = [];
    EEG = eeg_checkset(EEG);
    EEG = pop_select(EEG, 'trial', setdiff([1:EEG.trials], [1:4:EEG.trials]));
end;    

% this section adds ICA algorithms to the current path
% ----------------------------------------------------
addpath('/home/arno/matlab/allicas/02_FastICA_21');
addpath('/home/arno/matlab/allicas/05_JnS-1.2');
addpath('/home/arno/matlab/allicas/06_Pearson_ICA');
addpath('/home/arno/matlab/allicas/07_EGLD_ICA');
addpath('/home/arno/matlab/allicas/08_eeApack');
addpath('/home/arno/matlab/allicas/11_TFBSSpack');
addpath('/home/arno/matlab/allicas/14_icalab');
addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaMF');
addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaML');
addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaMS');
addpath('/data/common/matlab/eeglab/functions/sigprocfunc/');

% algorithms and parameters
% -------------------------
clear allalgs;
allalgs(1    ).algo = 'runica';        allalgs(end).speed = 517;   allalgs(end).name = 'Infomax';               allalgs(end).options = {};
allalgs(end+1).algo = 'runica';        allalgs(end).speed = 517;   allalgs(end).name = 'Ext. Infomax';          allalgs(end).options = { 'extended' 1 };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'FastICA';                      allalgs(end).options = {};
allalgs(end+1).algo = ''; %'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'fastica';    allalgs(end).options = { 'approach' 'symm' 'stabilization' 'on' };
allalgs(end+1).algo = ''; %'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'fastica';    allalgs(end).options = { 'approach' 'symm' 'stabilization' 'on' 'g' 'tanh' };
allalgs(end+1).algo = ''; % 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'fastica';    allalgs(end).options = { 'approach' 'symm' 'stabilization' 'on' 'g' 'gauss' };
allalgs(end+1).algo = ''; % 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'fastica';    allalgs(end).options = { 'approach' 'symm' 'stabilization' 'on' 'g' 'skew' };
allalgs(end+1).algo = 'jader';         allalgs(end).speed = 10668; allalgs(end).name = 'JADE';                  allalgs(end).options = {};
allalgs(end+1).algo = 'jadeop';        allalgs(end).speed = 9829;  allalgs(end).name = 'JADE opt.';             allalgs(end).options = {};
allalgs(end+1).algo = 'jade_td_p';     allalgs(end).speed = 765;   allalgs(end).name = 'JADE-TD';               allalgs(end).options = {};
allalgs(end+1).algo = 'MatlabshibbsR'; allalgs(end).speed = 4306;  allalgs(end).name = 'SHIBBS';                allalgs(end).options = {};
allalgs(end+1).algo = 'tica';          allalgs(end).speed = 615;   allalgs(end).name = 'TICA';                  allalgs(end).options = {};
allalgs(end+1).algo = 'erica';         allalgs(end).speed = 1053;  allalgs(end).name = 'ERICA';                 allalgs(end).options = {};
allalgs(end+1).algo = 'simbec';        allalgs(end).speed = 2611;  allalgs(end).name = 'simbec';                allalgs(end).options = {}; % Bad log-like, dipoles
allalgs(end+1).algo = 'unica';         allalgs(end).speed = 2509;  allalgs(end).name = 'unica';                 allalgs(end).options = {}; % Bad dipole curve
allalgs(end+1).algo = 'amuse';         allalgs(end).speed = 5;     allalgs(end).name = 'AMUSE';                 allalgs(end).options = {};
allalgs(end+1).algo = 'fobi';          allalgs(end).speed = 6;     allalgs(end).name = 'FOBI';                  allalgs(end).options = {};
allalgs(end+1).algo = 'evd';           allalgs(end).speed = 18;    allalgs(end).name = 'EVD';                   allalgs(end).options = {};
allalgs(end+1).algo = 'evd24';         allalgs(end).speed = 35;    allalgs(end).name = 'EVD24';                 allalgs(end).options = {};
allalgs(end+1).algo = 'sons';          allalgs(end).speed = 63717; allalgs(end).name = 'SONS';                  allalgs(end).options = {};
allalgs(end+1).algo = 'sobi';          allalgs(end).speed = 53;    allalgs(end).name = 'SOBI';                  allalgs(end).options = {};
allalgs(end+1).algo = 'ng_ol';         allalgs(end).speed = 955;   allalgs(end).name = 'nat. grad.';            allalgs(end).options = {}; % Do not know what it is
allalgs(end+1).algo = 'acsobiro';      allalgs(end).speed = 488;   allalgs(end).name = 'SOBIRO';                allalgs(end).options = {};
allalgs(end+1).algo = 'pearson_ica';   allalgs(end).speed = 4800;  allalgs(end).name = 'Pearson';               allalgs(end).options = {};
allalgs(end+1).algo = 'eeA';           allalgs(end).speed = 1956;  allalgs(end).name = 'eeA';                   allalgs(end).options = {};
allalgs(end+1).algo = 'icaML';         allalgs(end).speed = 8151;  allalgs(end).name = 'icaML';                 allalgs(end).options = {};
allalgs(end+1).algo = 'icaMS';         allalgs(end).speed = 138;   allalgs(end).name = 'icaMS';                 allalgs(end).options = {};
allalgs(end+1).algo = 'egld_ica';      allalgs(end).speed = 84883; allalgs(end).name = 'Egld';                  allalgs(end).options = {}; % never computed to completion

countalgs= length(allalgs);
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    allalgs(end).options = { 'approach' 'symm' }; % the one to show
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'g' 'tanh' 'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'g' 'gauss' 'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'g' 'skew' 'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'g' 'tanh'  'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'g' 'gauss' 'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'g' 'skew'  'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'epsilon' 0.00001  'displayMode' 'off' };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'epsilon' 0.000001 'displayMode' 'off'  };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'epsilon' 0.0000001 'displayMode' 'off'  };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'epsilon' 0.00000001 'displayMode' 'off'  };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'maxFineTune' 1000 'displayMode' 'off'  };
allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = '';    
allalgs(end).options = { 'approach' 'symm' 'maxNumIterations' 10000 'displayMode' 'off'  };
allalgs(end+1).algo = 'promax';        allalgs(end).speed = 0;     allalgs(end).name = 'Promax';    allalgs(end).options = {  };
allalgs(end+1).algo = 'pca';           allalgs(end).speed = 0;     allalgs(end).name = 'PCA';       allalgs(end).options = {  };
allalgs(end+1).algo = 'sphere';        allalgs(end).speed = 0;     allalgs(end).name = 'Sphering';  allalgs(end).options = {  };
allalgs(end+1).algo = 'formica';       allalgs(end).speed = 0;     allalgs(end).name = 'Amica';     allalgs(end).options = {  };
allalgs(end+1).algo = 'binica';        allalgs(end).speed = 0;     allalgs(end).name = 'binica';    allalgs(end).options = {  };
allalgs(end+1).algo = 'binica';        allalgs(end).speed = 0;     allalgs(end).name = 'binica ext.';    allalgs(end).options = { 'extended' 1 };

% ADD YOUR ALGORITHM HERE **************************************
%allalgs(48).algo = 'matlab_function';   % name of Matlab function
%allalgs(48).speed = 0;                  % speed of algorithm (simply enter 0)
%allalgs(48).name = 'This my ICA';       % name of Algorithm (for figures etc...)
%allalgs(48).options = { };              % options to give to the Matlab function (in addition to the data)

tmpind = cellfun('isempty', { allalgs.algo });
for ind = find(tmpind), allalgs(ind).name  = ''; end;

if exist('DATASET') ~= 1
    return;
end;
if exist('ALGONUM') ~= 1
    return;
end;

% tfbss = remove: did not converge with this data
% acrsobibpf was removed because it returned an empty matrix

filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
if ~exist(filename)
    if strcmpi(allalgs(ALGONUM).algo, 'sphere')
        tic; W = 2.0*inv(sqrtm(double(cov(EEG.data(:,:)')))); timeelapsed = toc;
    elseif strcmpi(allalgs(ALGONUM).algo, 'promax')
        tic; W = promax(EEG.data(:,:)); timeelapsed = toc;
    elseif strcmpi(allalgs(ALGONUM).algo, 'pca')
        tic; [pc Winv] = runpca(EEG.data(:,:)); timeelapsed = toc;
        W = inv(Winv);
    else
        tic; EEG = pop_runica(EEG, allalgs(ALGONUM).algo, allalgs(ALGONUM).options{:} ); timeelapsed = toc;
        W = EEG.icaweights*EEG.icasphere;
    end;
    eval( sprintf('save -mat %s W timeelapsed',filename));
else 
    load(filename);
end;

EEG.icawinv    = inv(W);
EEG.icaweights = W; 
EEG.icasphere  = eye(71); 

% use old dipfit plugin
tmpp = fileparts(which('dipplot'));
rmpath(tmpp);
tmpp = fileparts(which('processdat'));
addpath(fullfile(tmpp, 'dipfit1.02'));
addpath(fullfile(tmpp, 'dipfit1.02', 'copyprivate'));

EEG = eeg_multifit(EEG);
allrv = cell2mat( { EEG.dipfit.model.rv } );
for indexdip = 1:length(EEG.dipfit.model)
    allposxyz(:,indexdip) = EEG.dipfit.model(indexdip).posxyz(1,:)';
    allmomxyz(:,indexdip) = EEG.dipfit.model(indexdip).momxyz(1,:)';
    if all(allposxyz(:,indexdip) == 0), allposxyz(:,indexdip) = [NaN NaN NaN]; end;
    if all(allmomxyz(:,indexdip) == 0), allmomxyz(:,indexdip) = [NaN NaN NaN]; end;
end;
save( filename, '-mat', 'W', 'timeelapsed', 'allrv', 'allposxyz', 'allmomxyz');
return;

% launch on cluster
% -----------------
ALLALGONUM = length(allalgs):-1:1
%ALLALGONUM = [21 23 3 4 5 6 7 29:41 ]
ALLALGONUM = [4 24 26 27]
ALLALGONUM = [46 47]
fid = fopen('runallalgs', 'w');
for ALGONUM = ALLALGONUM
    for DATASET = 1:14
        if ~isempty(allalgs(ALGONUM).name)
            fprintf(fid, ['echo "cd matlab/mica/; matlab -nodesktop -nosplash -r ''' ...
                     'ALGONUM = ' int2str(ALGONUM) '; DATASET = ' int2str(DATASET) '; processdat; quit;''" | qsub -N ' ...
                     'mica' allalgs(ALGONUM).algo '_' int2str(ALGONUM) '_' int2str(DATASET) ' -o /home/arno/matlab/mica/qsub\n'], i, j);
        end;
    end;
end;
fclose(fid);
!chmod +x runallalgs

% test if file is present
% -----------------------
for DATASET = 1:14
    for ALGONUM = length(allalgs):-1:1
        if ~isempty(allalgs(ALGONUM).name)
            filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
            if ~exist(filename)
            fprintf(['echo "cd matlab/mica/; matlab -nodesktop -nosplash -r ''' ...
                     'ALGONUM = ' int2str(ALGONUM) '; DATASET = ' int2str(DATASET) '; processdat; quit;''" | qsub -N ' ...
                     'mica' allalgs(ALGONUM).algo '_' int2str(ALGONUM) '_' int2str(DATASET) ' -o /home/arno/matlab/mica/qsub -p -3\n'], i, j);
            end;
        end;
    end;
end;

% -----------------
% -----------------
% see plotresults for plotting results
% ------------------------------------



% plot each algorithm independanlty
% ---------------------------------
figure;
for ALGONUM = 1:length(allalgs)
    subplot(5,9, ALGONUM);
    plot(rvalgo(ALGONUM,:));
    title([ allalgs(ALGONUM).algo ' ' int2str(ALGONUM) ]);
end;

% plot maps for one algorithm
% ---------------------------
ALGONUM = 15;
DATASET = 1;
filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
load('-mat', filename);
figure; 
tmpinv = pinv(W);
for index = 1:20
    subplot(5,4,index);
    topoplot(tmpinv(:,index), EEG2.chanlocs);
    %topoplot(tmpinv(:,index), 'chan71.loc');
end;

% 5 curves
% --------
rvalgo = zeros(5, 71);
ALGONUM = 1;
datrange = 6:10;
for DATASET = datrange
    try
        filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
        load('-mat', filename);
        rvalgo(DATASET,:) = allrv;
    catch, 
        disp( [ 'Skipping ' int2str(DATASET) ' algo ' allalgs(ALGONUM).algo ] );
    end;
    rvalgo(DATASET,:) = sort(rvalgo(DATASET,:));
end;
figure; 
plot(rvalgo');
legend({'1' '2' '3' '4' '5'})
