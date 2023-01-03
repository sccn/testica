addpath './eeglab_current/eeglab2022.1'
ALGONUM = 48;
time = [];
remainingPmi = [];
dipolarity = [];
maxIter = 1000000;
%rules = [1e-8, 1e-9, 1e-10, 1e-11];
rules = [1e-1,1e-2,1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8];
nchans = 71;





for dat=1:14
   EEG = geticadata(dat, 'Picard');
   h0(:,dat) = getent2(reshape(EEG.data,nchans,EEG.pnts*EEG.trials));
end


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
        allalgs(end+1).algo = 'picard';        allalgs(end).speed = 0;     allalgs(end).name = 'Picard';    allalgs(end).options = { 'maxiter', maxIter, 'mode', 'standard' };
        allalgs(end+1).algo = 'picard';        allalgs(end).speed = 0;     allalgs(end).name = 'PicardO';    allalgs(end).options = { 'maxiter', maxIter, 'mode', 'ortho' };
        tmpind = cellfun('isempty', { allalgs.algo });
        for ind = find(tmpind), allalgs(ind).name  = ''; end;


PMI = zeros(nchans,nchans,length(rules),14);
mir = zeros(length(rules),14);
h = zeros(71,length(rules),14);
totalDipole = zeros(length(rules),14);
completeRV = zeros(71,length(rules),14);
for stoppingRuleIdx=1:length(rules)
    stoppingRule = rules(stoppingRuleIdx);
    sRDipolarity = [];
    sRTotalPMI = [];
    parfor d=1:14
    %for d = 1
        % have to define 2 values% have to define 2 values
        % DATASET: sample dataset to use
        % ALGONUM: algorythm number
        DATASET = d;
        if ~exist('pop_loadset.m','file')
            eegladipolarity = [];b; close;
        end;
        %clear ALLEEG EEG CURRENTSET ALLCOM LASTCOM;
        
        parentpath = fileparts(which('/processdat'));
        parentpath = fullfile(parentpath, 'datasets');
        if true
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
        addpath('./eeglab_current/eeglab2021.1/functions/sigprocfunc/')
        
        % algorithms and parameters
        % -------------------------
       
        
      
        
        % tfbss = remove: did not converge with this data
        % acrsobibpf was removed because it returned an empty matrix
        
        %filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
        tic; EEG = pop_runica(EEG, allalgs(ALGONUM).algo, allalgs(ALGONUM).options{:}, 'tol', stoppingRule); timeelapsed = toc;
        W = EEG.icaweights*EEG.icasphere;


        s = W * reshape(EEG.data,nchans,EEG.pnts*EEG.trials); 
        [MI,T,Hu,hu,Hx] = get_mi(s);
        %get PMI
        PMI(:,:,stoppingRuleIdx,d) = MI;
        %get mir
        h(:,stoppingRuleIdx,d) = getent2(s);
        mir(stoppingRuleIdx,d) = sum(h0(:,d)) - sum(h(:,stoppingRuleIdx,d)) + sum(log(abs(eig(W))));

        EEG.icawinv    = inv(W);
        EEG.icaweights = W; 
        EEG.icasphere  = eye(71); 
        
        % use old dipfit plugin
        tmpp = fileparts(which('/dipplot'));
        rmpath(tmpp);
        tmpp = fileparts(which('/processdat'));
        addpath(fullfile(tmpp, 'dipfit1.02'));
        addpath(fullfile(tmpp, 'dipfit1.02', 'copyprivate'));
        
        EEG = eeg_multifit(EEG);
        allrv = cell2mat( { EEG.dipfit.model.rv } );

%         for indexdip = 1:length(EEG.dipfit.model)
%             allposxyz(:,indexdip) = EEG.dipfit.model(indexdip).posxyz(1,:)';
%             allmomxyz(:,indexdip) = EEG.dipfit.model(indexdip).momxyz(1,:)';
%             if all(allposxyz(:,indexdip) == 0), allposxyz(:,indexdip) = [NaN NaN NaN]; end;
%             if all(allmomxyz(:,indexdip) == 0), allmomxyz(:,indexdip) = [NaN NaN NaN]; end;
%         end;

        %store allrv
        completeRV(:,stoppingRuleIdx,d) = allrv;
        %sRDipolarity(end+1) = allrv;z
        %sRDipolarity(end+1) = length(find(allrv < 0.05));
        totalDipole(stoppingRuleIdx, d) = length(find(allrv < 0.05));
        %save( filename, '-mat', 'W', 'timeelapsed', 'allrv', 'allposxyz', 'allmomxyz');
        %return;
        
    end
    %dipolarity(end+1) = mean(sRDipolarity);
    filename = strcat('dipolarityPMIMIRRulesSept26_',int2str(stoppingRuleIdx),'.mat');
    save( filename, '-mat', 'totalDipole','PMI','mir');
end
filename = 'dipolarityPMIMIRRulesSept26Final.mat';
save( filename, '-mat', 'totalDipole','PMI','mir');
%plot PMI
%load('dipolarityPMIRes.mat')
%calculate PMI for raw EEG
PMIraw = zeros(nchans,nchans,14);
totalPMIraw = zeros(14);

%get the raw PMI
for dat=1:14
   EEG = geticadata(dat, 'Picard');
   [MI,T,Hu,hu,Hx] = get_mi(reshape(EEG.data,nchans,EEG.pnts*EEG.trials));
   PMIraw(:,:,dat) = MI;
   totalPMIraw = sum(MI);
end
%dataset 10 is an outlier
goodSub = [1 2 3 4 5 6 7 8 9 11 12 13 14];
PMIrawslice = PMIraw(:,:,goodSub);
PMIslice = PMI(:,:,:,goodSub);
PMIperm = permute(PMIslice, [3 4 1 2]);
PMIperm = sum(PMIperm, [3,4]);
PMIrawperm = permute(PMIrawslice, [3 1 2]);
PMIrawperm = sum(PMIrawperm, [2,3]);
    %mir.mir = mir.PMI./repmat(mir2.PMIraw, [size(mir.PMI,1) 1])*100;
for i=1:length(rules) %number of algorithms
    PMIReduction(i,:) = PMIperm(i,:)./PMIrawperm'*100;
end

%totalPMI = squeeze(sum(PMIslice,[1,2]));
avgPMI = mean(PMIReduction,2);
stdPMI = std(PMIReduction,0,2);

figure('position', [10 10 900 900])
%xticklabels = {'Infomax';'Ext. Infomax';'FastICA';'Picard'}
asscendingRules = flip(rules);
asscendingAvgPMI = flip(avgPMI);
asscendingStdPMI = flip(stdPMI);
bh = bar(log(asscendingRules),asscendingAvgPMI);                
hold on
err = errorbar(log(asscendingRules),asscendingAvgPMI,asscendingStdPMI,'lineWidth',2);
err.Color = [0 0 0];    
err.LineStyle = 'none'; 
bh.FaceColor = 'flat';
bh.CData = [0 0 1; 0 0 1; 0 0 1; 1 0 0];
xticklabels = asscendingRules;
set(gca,'xtick',log(asscendingRules),'xticklabel',xticklabels)
xlabel('Stopping Rule')
ylabel('Remnant PMI')
print('PMI_rules.eps','-depsc')
hold off

figure('position', [10 10 900 900])
%xticklabels = {'Infomax';'Ext. Infomax';'FastICA';'Picard'}
mirSlice = mir([1 2 3 4 5 6],:);
avgMIR = mean(mirSlice,2);
asscendingRules = flip(rules([1 2 3 4 5 6]));
asscendingAvgMIR = flip(avgMIR);
%asscendingStdPMI = flip(stdPMI);
bh = bar(log(asscendingRules),asscendingAvgMIR);                
hold on
%err = errorbar(log(asscendingRules),asscendingAvgPMI,asscendingStdPMI,'lineWidth',2);
%err.Color = [0 0 0];    
%err.LineStyle = 'none'; 
%bh.FaceColor = 'flat';
%bh.CData = [0 0 1; 0 0 1; 0 0 1; 1 0 0];
xticklabels = asscendingRules;
set(gca,'xtick',log(asscendingRules),'xticklabel',xticklabels)
xlabel('Stopping Rule')
ylabel('Mutual Information Reduction')
print('MIR_rules.eps','-depsc')
hold off

%plot dipolarity
dipoleSlice = totalDipole([1 2 3 4 5 6],[1 2 3 4 5 6 7 8 9 11 12 13 14])
avgDipole = mean(dipoleSlice,2)
figure('position', [10 10 900 900])
%xticklabels = {'Infomax';'Ext. Infomax';'FastICA';'Picard'}
asscendingRules = flip(rules([1 2 3 4 5 6]));
asscendingDipolarity = flip(avgDipole);
bh = bar(log(asscendingRules),asscendingDipolarity);                
hold on
xticklabels = asscendingRules;
set(gca,'xtick',log(asscendingRules),'xticklabel',xticklabels)
xlabel('Stopping Rule')
ylabel('Percent near-dipolar (r.v. < 5%) components')
print('dipolarity_rules.eps','-depsc')
hold off
