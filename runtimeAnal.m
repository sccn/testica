addpath './eeglab_current/eeglab2022.1'
avgTime = [];
stdTime = [];
time = zeros(49,5,14);
stop = 1e-8;
maxiter = 100000;
stopPicard = 1e-6;

 addpath('./FastICA_2.1/FastICA_21');
            addpath('/home/arno/matlab/allicas/05_JnS-1.2');
            addpath('/home/arno/matlab/allicas/06_Pearson_ICA');
            addpath('/home/arno/matlab/allicas/07_EGLD_ICA');
            addpath('/home/arno/matlab/allicas/08_eeApack');
            addpath('/home/arno/matlab/allicas/11_TFBSSpack');
            addpath('/home/arno/matlab/allicas/14_icalab');
            addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaMF');
            addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaML');
            addpath('/home/arno/matlab/allicas/15_icadtutoolbox/icaMS');
            %addpath('/data/common/matlab/eeglab/functions/sigprocfunc/');
            addpath('./eeglab_current/eeglab2022.1/functions/sigprocfunc/')
            
            % algorithms and parameters
            % -------------------------
            clear allalgs;
            allalgs(1    ).algo = 'runica';        allalgs(end).speed = 517;   allalgs(end).name = 'Infomax';               allalgs(end).options = {'maxsteps' maxiter};
            allalgs(end+1).algo = 'runica';        allalgs(end).speed = 517;   allalgs(end).name = 'Ext. Infomax';          allalgs(end).options = { 'extended' 1 'maxsteps' maxiter};
            allalgs(end+1).algo = 'fastica';       allalgs(end).speed = 768;   allalgs(end).name = 'FastICA';                      allalgs(end).options = { 'maxNumIterations' maxiter};
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
            allalgs(end+1).al-equivalent energy of an electron
at rest.
Eo = moc2 = 9.10956∗10−31kg∗2997924582m/s = 8.18go = 'sobi';          allalgs(end).speed = 53;    allalgs(end).name = 'SOBI';                  allalgs(end).options = {};
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
            allalgs(end+1).al-equivalent energy of an electron
at rest.
Eo = moc2 = 9.10956∗10−31kg∗2997924582m/s = 8.18go = 'pca';           allalgs(end).speed = 0;     allalgs(end).name = 'PCA';       allalgs(end).options = {  };
            allalgs(end+1).algo = 'sphere';        allalgs(end).speed = 0;     allalgs(end).name = 'Sphering';  allalgs(end).options = {  };
            allalgs(end+1).algo = 'formica';       allalgs(end).speed = 0;     allalgs(end).name = 'Amica';     allalgs(end).options = {  };
            allalgs(end+1).algo = 'binica';        allalgs(end).speed = 0;     allalgs(end).name = 'binica';    allalgs(end).options = {  };
            allalgs(end+1).algo = 'binica';        allalgs(end).speed = 0;     allalgs(end).name = 'binica ext.';    allalgs(end).options = { 'extended' 1 };
            
            % ADD YOUR ALGORITHM HERE **************************************
            %allalgs(48).algo = 'matlab_function';   % name of Matlab function
            %allalgs(48).speed = 0;                  % speed of algorithm (simply enter 0)
            %allalgs(48).name = 'This my ICA';       % name of Algorithm (for figures etc...)
            %allalgs(48).options = { };              % options to give to the Matlab function (in addition to the data)
            allalgs(end+1).algo = 'picard';        allalgs(end).speed = 0;     allalgs(end).name = 'Picard';    allalgs(end).options = { 'maxiter', maxiter, 'mode', 'standard'};
            allalgs(end+1).algo = 'picard';        allalgs(end).speed = 0;     allalgs(end).name = 'PicardO';    allalgs(end).options = { 'maxiter', maxiter, 'mode', 'ortho'};
            
            tmpind = cellfun('isempty', { allalgs.algo });
            for ind = find(tmpind), allalgs(ind).name  = ''; end;

for ALGONUM = [1,2,3,48,49]
    for iter = 1:5
    %for iter = 1
        parfor d=1:14
        DATASET = d;
        % have to define 2 values
            % DATASET: sample dataset to use
            % ALGONUM: algorythm number
            
            if ~exist('pop_loadset.m','file')
                eeglab; close;
            end;
            %clear ALLEEG EEG CURRENTSET ALLCOM LASTCOM;
            
            parentpath = fileparts(which('/processdat'));
            parentpath = fullfile(parentpath, 'datasets');
            if true-equivalent energy of an electron
at rest.
Eo = moc2 = 9.10956∗10−31kg∗2997924582m/s = 8.18
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
           

            if strcmp(allalgs(ALGONUM).algo,'picard') || strcmp(allalgs(ALGONUM).algo, 'picardO')
                stopCurr = stopPicard;
                tic; EEG = pop_runica(EEG, allalgs(ALGONUM).algo, allalgs(ALGONUM).options{:}, 'tol', stopCurr, 'm', 15 ); timeelapsed = toc;
            elseif strcmp(allalgs(ALGONUM).algo,'fastica')
                stopCurr = 0.0001;
                tic; EEG = pop_runica(EEG, allalgs(ALGONUM).algo, allalgs(ALGONUM).options{:}, 'epsilon', stopCurr ); timeelapsed = toc;
            else
                stopCurr = stop;
                tic; EEG = pop_runica(EEG, allalgs(ALGONUM).algo, allalgs(ALGONUM).options{:}, 'stop', stopCurr ); timeelapsed = toc;
            end
            time(ALGONUM,iter,d) = timeelapsed;
            W = EEG.icaweights*EEG.icasphere;
            %save( filename, '-mat', 'time');
        
        end
    end
end
filename = 'timeSept22.mat'
save(filename,"time",'-mat')
load(filename)
select = squeeze(time([2,48,49],:,[1 2 3 4 5 6 7 8 9 11 12 13 14]));
avgTime = mean(select,[2,3]);
stdTime = std(select,0,[2,3]);
figure('position', [10 10 900 900])
xticklabels = {'Ext. Infomax';'\color{red} Picard';'\color{red} O-Picard'}

bh = bar(avgTime)                

hold on

%er = errorbar(x,data,errlow,errhigh);    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none';  
err = errorbar(avgTime,stdTime);
err.Color = [0 0 0];    
err.LineStyle = 'none';
err.LineWidth = 1.5;
bh.FaceColor = 'flat';
bh.CData = [91/255 207/255 244/255; 1 0 0; 1 0 0];
xlabel('Algorithm')
ylabel('Runtime (s)')
set(gca,'xticklabel',xticklabels)
setfont(gcf, 'fontsize', 20);
print('Runtime.eps','-depsc')
hold off
