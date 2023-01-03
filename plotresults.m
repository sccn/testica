% This script plots results for all ICA algorithm
% Because some ICA algorithm are computed and not plotted
% you must add the name of your algorithm to the list
% Arnaud Delorme - Aug 2011

clear;
newalgs = {  'Amica' 'Ext. Infomax' 'Pearson' 'Infomax' ...
                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
                 'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
                'AMUSE'  'Picard' 'Picard-O' 'PCA'}; % 'icaML' 'promax' 'unica'    'simbec' 
newalgsColor = {  'Amica' 'Ext. Infomax' 'Pearson' 'Infomax' ...
                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
                 'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
                'AMUSE'  '\color{red} Picard' '\color{red} Picard-O' 'PCA'};

plotmir = 1; % plot MIR versus PMI

% -----------------
% end of parameters
% -----------------

if exist('datasetrange')
    if ~exist('ALLEEG')
        ALLEEG = [];
    end;
    set(0, 'userdata', { datasetrange ALLEEG } );
end;
processdat; % load allalgs variable

% reload temporary saved variables
userdata = get(0, 'userdata');
if isempty(userdata)
    datasetrange = [1:9, 11:14]; %1:14;
else
    datasetrange = userdata{1};
    ALLEEG       = userdata{2};
end;

nbdat        = length(datasetrange);
inc          = 1/nbdat;

% import results and plot
% -----------------------
allres.rvalgo     = zeros(length(allalgs), 71*nbdat);
allres.allpos     = zeros(length(allalgs), 3,  71*nbdat);
allres.allmom     = zeros(length(allalgs), 3,  71*nbdat);
allres.winv       = zeros(length(allalgs), 71, 71*nbdat);
allres.algtime    = zeros(1,length(allalgs));
clear allres;

countplot  = 1;
for ind = length(newalgs):-1:1
    addtoplot    = 0;
    countdat     = 1;
    ALGONUM      = strmatch(newalgs{ind}, { allalgs.name }, 'exact');
    
    if ~isempty(newalgs{ind}) & ~isempty(ALGONUM)
        for DATASET = datasetrange
            try
                clear allrv
                filename = sprintf('icadecompositions/ica%d_72_%2.2d_%s.mat',DATASET, ALGONUM, allalgs(ALGONUM).algo);
                load('-mat', filename);
                allres(countplot).rvalgo (:, countdat)     = allrv;
                allres(countplot).allpos (:, :, countdat)  = allposxyz;
                allres(countplot).allmom (:, :, countdat)  = allmomxyz;

                allres(countplot).winv(:, :, countdat)     = pinv(W);
                alltime   (countdat)                       = timeelapsed;
                allbelow10(countdat)                       = length(find(allrv < 0.1));
                addtoplot = addtoplot+1;

                % Plotting result of a specific algo and dataset
                % ----------------------------------------------
                %if ALGONUM == 2 && DATASET == 2
                %    figure; 
                %    allallres.winv = pinv(W);
                %    for index = 1:20
                %        subplot(5,4,index);
                %        topoplot(allallres.winv(:,index,1), 'chan71.loc');
                %        title(['ext. runica ' num2str(allrv(index),3) ]);
                %    end;
                %    return;
                %end;

            catch, 
                disp(lasterr);
                disp( [ 'Skipping ' int2str(DATASET) ' algo ' allalgs(ALGONUM).algo ' ' lasterr ] );
            end;
            countdat = countdat + 1;
        end;
    end;
    if addtoplot == length(datasetrange)
        %if strcmp(allalgs(ALGONUM).name, 'Picard') | strcmp(allalgs(ALGONUM).name, 'Picard-O')
        %    allres(countplot).legend   = strcat('\color{red} ',allalgs(ALGONUM).name);
        %else
            allres(countplot).legend   = allalgs(ALGONUM).name;
        %end
        allres(countplot).alginfo  = allalgs(ALGONUM);
        allres(countplot).algtime  = mean(alltime);
        fprintf( 'Algo %s: %s\n', allalgs(ALGONUM).algo, int2str(allbelow10) );

        allres(countplot).std_time      = std (alltime);
        allres(countplot).stdrvbelow10p = std (allbelow10);
        countplot = countplot+1;
    else
        newalgs(ind) = [];
        try,
            allres(countplot) = [];
        catch, end;
    end;
end;

% sort algorithms by performance (number of components for residual variance < 10%)
% ---------------------------------------------------------------------------------
for index = 1:length(allres)
    allres(index).allrvsorted = sort(allres(index).rvalgo(:));
    allres(index).rvbelow10   = length(find(allres(index).rvalgo(:) < 0.1));
    allres(index).rvbelow4    = length(find(allres(index).rvalgo(:) < 0.04));
    allres(index).rvbelow5    = length(find(allres(index).rvalgo(:) < 0.05));
    allres(index).rvbelow6    = length(find(allres(index).rvalgo(:) < 0.06));
    allres(index).rvbelow15   = length(find(allres(index).rvalgo(:) < 0.15));
    allres(index).rvbelow75   = length(find(allres(index).rvalgo(:) < 0.075));
    allres(index).rvbelow10p  = allres(index).rvbelow10/length(datasetrange);
    
    for s = 1:size(allres(index).rvalgo,2)
        allres(index).rvbelow4subj(s) = length(find(allres(index).rvalgo(:,s) < 0.04));
        allres(index).rvbelow5subj(s) = length(find(allres(index).rvalgo(:,s) < 0.05));
        allres(index).rvbelow6subj(s) = length(find(allres(index).rvalgo(:,s) < 0.06));
        allres(index).rvbelow10subj(s) = length(find(allres(index).rvalgo(:,s) < 0.1));
    end;
    
    % get rv at all threshold by 0.01 increments
    for percent = 1:100
        allres(index).rvbelowp(percent) = length(find(allres(index).rvalgo(:) < percent/100));
    end;
end;

% reorder algos
% -------------
if 0
    [tmp so]     = sort([ allres.rvbelow10 ]);
    [tmp sorev]  = sort(so(end:-1:1)); 
    allres       = allres(so(end:-1:1));
else
    allres = allres(end:-1:1);
end;

% define 13 colors
% ----------------
allcolors = hsv(12);
allcolors(1 ,:) = [1.0000    0.6855    0.1109];
allcolors(2 ,:) = [0.5 0.5 0.5];
allcolors(4 ,:) = [0   0.6 0];
allcolors(6 ,:) = [0.8 0.8 0];
allcolors(9 ,:) = [0   0   0];
allcolors(13,:) = [0.7 0   0.7];
allcolors(14, :) = [0.7 1 0.7];
allcolors(3,:)  = [];
%allcolors = hsv(22);
%figure; imagesc(reshape(allcolors, size(allcolors,1), 1, 3));

% add EEG in the background
% -------------------------
clear totallrv
count = 1;
for index = datasetrange
    filename2 = sprintf('randdip/randomdip_%d_data.mat',index);
    load('-mat', filename2);
    totallrv(count,:) = cumout(:)';
    count = count+1;
end;
totallrv = totallrv';
tmpallrv = reshape(totallrv, 71, 10, size(totallrv,2));
totallrv = sort(totallrv,1);
%readploteeg;

% plot result of the algorithms
% -----------------------------
%first let's make a sorted version of allres so we can plot algorithms in 
%descending order of perecent N.D. components
T = struct2table(allres);
sortedT = sortrows(T, 'rvbelow10','descend');
allresSorted = table2struct(sortedT);

figure('position', [10   10   1000   1000]);
hh = semilogy(allresSorted(1).allrvsorted, 'w'); hold on;
delete(hh);
for index = 1:length(allresSorted)
    %tmpind = find( allres.rvalgo(:,index) );
    %tmpdat = allres.rvalgo(tmpind,index);
    %tmpind = ~isnan(tmpdat);
    %tmpdat = tmpdat(tmpind);
    %abscicia = linspace(1,72, length(tmpdat));
    %if length(tmpdat) ~= 0
    %   h = semilogy(abscicia, tmpdat); hold on;
    h = semilogy(linspace(1,100,length(allresSorted(index).allrvsorted)), allresSorted(index).allrvsorted); hold on;
    %h = plot([1:inc:71.8], log10(allres.rvalgo(:,index))); hold on;
    set(h, 'linewidth', 2, 'color', allcolors(mod(index-1,13)+1,:));
    if contains(allresSorted(index).legend, 'Picard-O') %we'll pay extra attention to Picard and make sure to highlight the results
        set(h, 'linestyle', '-', 'LineWidth',8,'Color',[1 .31 0]);
    elseif contains(allresSorted(index).legend, 'Picard') %we'll pay extra attention to Picard and make sure to highlight the results
        set(h, 'linestyle', ':', 'LineWidth',8,'Color',[0 .69 1]);
    elseif strcmp(allresSorted(index).legend, 'PCA') %we'll pay extra attention to Picard and make sure to highlight the results
        set(h, 'linestyle', ':', 'LineWidth',8, 'Color',[1 0 1]);
    else
        if index>12, set(h, 'linestyle', '--'); end;
    end
end;
%figure; plot(sqrt([1:71*5]), log10(allres.rvalgo'));
%figure; semilogx(log10(allres.rvalgo'));
hh = semilogy(linspace(1,100, 710), nan_mean(nan_mean(totallrv,2),2), 'k-.', 'linewidth', 4); hold on;
lh = legend({ allresSorted.legend 'raw EEG' });
view([90 270]);
xlabel('Percent of ICA components');
ylabel('Dipole model residual variance (%)');
set(gca, 'ylim', [0.004, 1]);
set(gca, 'xlim', [0 100]);
set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
set(gcf, 'paperpositionmode', 'auto');
setfont(gcf, 'fontsize', 20);
set(lh, 'position', [0.1598 0.2119 0.1670 0.7141]);
%eval(sprintf('print -depsc icacomp%d-%d.eps', datasetrange(1), datasetrange(end)));
%print -djpeg icacomp.jpg

saveas(gcf,'rv_vs_percentICA.png')

% % plot PMI versus MIR
% -------------------
clear allmir allpmi;
mir = load( '-mat', 'pmi_new.mat');
mir2 = load('-mat', 'pmi_new_raw.mat');
mir.PMI = permute(mir.PMI, [3 4 1 2]);
mir.PMI = sum(mir.PMI, [3,4]);
mir2.PMIraw = permute(mir2.PMIraw, [3 1 2]);
mir2.PMIraw = sum(mir2.PMIraw, [2,3]);
for i=1:20 %number of algorithms
   mir.mir(i,:) = mir.PMI(i,:)./mir2.PMIraw'*100;
end

pmi = mir;
mir = load('-mat', 'mir_new.mat');
mir.mir = mir.mir*1.4427*250;
figure;
for dat = 1:length(datasetrange)
    for algo=1:length(newalgs)
        inda = strmatch(lower(newalgs{algo}), lower(mir.algorithms), 'exact');
        indb = strmatch(lower(newalgs{algo}), lower(pmi.algorithms), 'exact');
        allmir(algo) = mean(mir.mir(inda,datasetrange(dat))-mir.mir(end-1,datasetrange(dat)));
        allpmi(algo) = mean(pmi.PMI(indb,datasetrange(dat))-pmi.PMI(end,datasetrange(dat)));
    end;
    %allresmir(dat) = (allmir(1)/allmir(end-1)-1)*100; % (MIR_AMICA/MIR_PCA-1)*100
    %allrespmi(dat) = (allpmi(1)/allpmi(end-1)-1)*100;    
    plot(allmir, allpmi, '.', 'color', allcolors(dat,:));
    [ypred alpha10 r10 slope] = myregress(allmir, allpmi);
    hold on; plot(allmir, ypred, 'k:', 'linewidth', 2, 'color', allcolors(dat,:));
end;
xlabel('MIR');
ylabel('PMI');

%% plot MI or PMI versus number of dipoles
% --------------------------------
%newalgs = { 'Amica'    'Infomax'    'Ext. Infomax'    'Pearson'    'SHIBBS'    'JADE'    'FastICA' 'TICA'    'JADE opt.'    'SOBI'    'JADE-TD'    'SOBIRO'  };
clear mir allpmi allmir alldip5 alldip10 alldip15 alldip75;
%pmi = load( '-mat', 'pmi_save.mat');
%pmi2 = load('-mat', 'pmi_save_data.mat');
%pmi.pmi = pmi.PMI./repmat(pmi2.PMI, [size(pmi.PMI,1) 1])*100;
pmiOld = load( '-mat', 'pmi_stat_save_oct2011.mat');
pmi2Old = load('-mat', 'pmi_stat_save_raw_oct2011.mat');
%pmi.pmi = pmi.PMI*1.4427*250*2845/1000;
%pmi.pmi = pmi.PMI./repmat(pmi2.PMIraw, [size(pmi.PMI,1) 1])*100; % in %
if plotmir % CHANGE HERE FOR PMI VERSUS MIR
    %plot MIR
    mir = load('-mat', 'mir_new.mat');
    mir.mir = mir.mir*1.4427*250; %Convert MIR to bits/sec
    %mir.std = std(mir.mir,0,2);
    %compute Standard Deviation
else 
    
    mir = load( '-mat', 'pmi_new.mat');
    mir2 = load('-mat', 'pmi_new_raw.mat');
    mir.PMI = permute(mir.PMI, [3 4 1 2]);
    mir.PMI = mean(mir.PMI, [3,4]);
    mir2.PMIraw = permute(mir2.PMIraw, [3 1 2]);
    mir2.PMIraw = mean(mir2.PMIraw, [2,3])';
    for i=1:length(mir.algorithms) %number of algorithms
        mir.mir(i,:) = (mir.PMI(i,:)./mir2.PMIraw)*100;
    end

end;
tmpname = cell(1,length(allres)); tmpname(:) = { 'name' };
algnames = cellfun(@getfield, { allres.alginfo }, tmpname, 'UniformOutput', false);
for algo=1:length(newalgs)
    inda   = strmatch(lower(newalgs{algo}), lower(mir.algorithms), 'exact');
    indb   = strmatch(lower(newalgs{algo}), lower(pmi.algorithms), 'exact');
    indres = strmatch(lower(newalgs{algo}), lower({ allres.legend }), 'exact');
    
    %allmir(:,algo) = (mir.mir(end,datasetrange)-mir.mir(inda,datasetrange))*1.4427; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits | 250 is the sampling rate
    allmir(:,algo) = mir.mir(inda,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    %%allpmi(:,algo) = pmi.pmi(indb,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    %%allpmi(:,algo) = pmi2.PMIraw(datasetrange)*1.4427*250*2845/1000; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    allres(indres).mir = allmir(:,algo);
    %%allres(indres).pmi = allpmi(:,algo);
    if isfield(mir, 'std')
        allmirstd(:,algo) = mir.std(inda,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    end;
    
    inda = strmatch(newalgs{algo}, algnames, 'exact');
    alldip4(algo) = allres(inda).rvbelow4/length(datasetrange)/71*100;
    alldip5(algo) = allres(inda).rvbelow5/length(datasetrange)/71*100;
    %alldip5dot5(algo) = allres(inda).rvbelow5dot5/length(datasetrange)/71*100;
    alldip6(algo) = allres(inda).rvbelow6/length(datasetrange)/71*100;
    alldip10(algo) = allres(inda).rvbelow10/length(datasetrange)/71*100;
    alldip15(algo) = allres(inda).rvbelow15/length(datasetrange)/71*100;
    alldip75(algo) = allres(inda).rvbelow75/length(datasetrange)/71*100;
    alldip4subj(algo,:)  = allres(inda).rvbelow4subj /71*100;
    alldip5subj(algo,:)  = allres(inda).rvbelow5subj /71*100;
    alldip6subj(algo,:)  = allres(inda).rvbelow6subj /71*100;
    alldip10subj(algo,:) = allres(inda).rvbelow10subj/71*100;
    alldipp(algo,:) = allres(inda).rvbelowp/length(datasetrange)/71*100';
end;
figure('position', [10 10 900 900]); hold on;
for i=1:length(alldip10)
    if strcmp(newalgs{i}, "Picard") | strcmp(newalgs{i}, "Picard-O")
        if plotmir == 1
            plot(mean(allmir(:,i),'omitnan')/1000, alldip5(i), '.', 'color', 'r', 'markersize', 18);
        else
            plot(mean(allmir(:,i),'omitnan'), alldip5(i), '.', 'color', 'r', 'markersize', 18);
        end
    else
        if plotmir == 1
            plot(mean(allmir(:,i),'omitnan')/1000, alldip5(i), '.', 'color', 'k', 'markersize', 18);
        else
            plot(mean(allmir(:,i),'omitnan'), alldip5(i), '.', 'color', 'k', 'markersize', 18);
        end
    end
    if plotmir == 1
        allmirstd(:,i) = std(allmir(:,i)/1000)*.2;
    else
        allmirstd(:,i) = std(allmir(:,i))*.2;
    end
    if 1%if isfield(mir, 'std')
        if plotmir == 1
            valx1 = (mean(allmir(:,i)/1000)-mean(allmirstd(:,i)));
            valx2 = (mean(allmir(:,i)/1000)+mean(allmirstd(:,i)));
        else
            valx1 = (mean(allmir(:,i))-mean(allmirstd(:,i)));
            valx2 = (mean(allmir(:,i))+mean(allmirstd(:,i)));
        end
        %perfectly so to center the ellipse we'll just center it at
        rvStd = std(alldip5subj(i,:))*.2;
        valy1 = alldip5(i)-rvStd;
        valy2 = alldip5(i)+rvStd;
       if strcmp(newalgs{i}, "Picard") | strcmp(newalgs{i}, "Picard-O")
        ellipse((valx2-valx1)/2, (valy2-valy1)/2, 0, (valx1+valx2)/2, (valy1+valy2)/2, 'r') %computes second axis using STD values
       else
        ellipse_dot((valx2-valx1)/2, (valy2-valy1)/2, 0, (valx1+valx2)/2, (valy1+valy2)/2, [.1 .1 .1]) %computes second axis using STD values
       end
       
    end;

    if plotmir
        if strcmp(newalgs{i}, "Picard") | strcmp(newalgs{i}, "Picard-O")
            text(mean(allmir(:,i))/1000+0.02, alldip5(i), newalgs{i},'color','r');
        else
            text(mean(allmir(:,i))/1000+0.02, alldip5(i), newalgs{i},'color','k');
        end
        
    else
        if strcmp(newalgs{i}, "Picard") | strcmp(newalgs{i}, "Picard-O")
            text(mean(allmir(:,i),'omitnan')+0.00002, alldip5(i), newalgs{i},'color','r');
        else
            text(mean(allmir(:,i),'omitnan')+0.00002, alldip5(i), newalgs{i},'color','k');
        end
        
    end;
end;

%% xlim([-0.05 2]*1.4427)
if plotmir
    [ypred alpha5 r5 slope] = myregress(mean(allmir)/1000, alldip5);
    fprintf('-----------------------------------\nRegression at 5%% threshold: r^2=%1.3f\tp=%e\n', r5, alpha5);
    hold on; plot(mean(allmir)/1000, ypred, 'k--', 'linewidth', 1.5);
    myStd = std(allmir)/1000
else
     [ypred alpha5 r5 slope] = myregress(mean(allmir), alldip5);
    fprintf('-----------------------------------\nRegression at 5%% threshold: r^2=%1.3f\tp=%e\n', r5, alpha5);
    hold on; plot(mean(allmir), ypred, 'k--', 'linewidth', 1.5);
    myStd = std(allmir)
end
if 0
    %for i=1:length(alldip10), plot(mean(allmir(:,i)), alldip75(i), '.', 'color', allcolors(i,:), 'markersize', 18); end;
    %[ypred alpha75 r75 slope] = myregress(mean(allmir), alldip75);
    %hold on; plot(mean(allmir), ypred, 'k', 'linewidth', 1.5);
    
    for i=1:length(alldip10), plot(mean(allmir(:,i)), alldip10(i), '.', 'color', allcolors(i,:), 'markersize', 18); end;
    myStd = std(allmir(:,i))
    [ypred alpha10 r10 slope] = myregress(mean(allmir), alldip10);
    hold on; plot(mean(allmir), ypred, 'k:', 'linewidth', 2);
end;
setfont(gcf, 'fontsize', 20);
if plotmir
    xlabel('Mutual information reduction (kbits/sec)');
else
    xlabel('Remnant Pairwise Mutual Information (%)', 'interpreter', 'none');
end;
ylabel('Percent near-dipolar (r.v. < 5%) components');

averageMir = mean(allmir)/1000
save('tableOut.mat','newalgs','averageMir','alldip5');

if plotmir
saveas(gcf,'MIR_vs_ND.png')
print('MIR_vs_ND','-depsc')
else
saveas(gcf,'PMI_vs_ND.png')
print('PMI_vs_ND.eps','-depsc')
end


%% plot MI or PMI versus number of dipoles
% --------------------------------
%newalgs = { 'Amica'    'Infomax'    'Ext. Infomax'    'Pearson'    'SHIBBS'    'JADE'    'FastICA' 'TICA'    'JADE opt.'    'SOBI'    'JADE-TD'    'SOBIRO'  };
clear mir allpmi allmir alldip5 alldip10 alldip15 alldip75;
%TODO: Edit to match your saved pmi data
pmi = load( '-mat', 'pmi_stat_save_oct2011.mat');
pmi2 = load('-mat', 'pmi_stat_save_raw_oct2011.mat');
pmi.pmi = pmi.PMI./repmat(pmi2.PMIraw, [size(pmi.PMI,1) 1])*100; % in %
if plotmir % CHANGE HERE FOR PMI VERSUS MIR
    mir = load('-mat', 'mir_new.mat');
    mir.mir = mir.mir*1.4427*250; %Convert MIR to bits/sec
else
    mir = load( '-mat', 'pmi_new.mat');
    mir2 = load('-mat', 'pmi_new_raw.mat');
    mir.PMI = permute(mir.PMI, [3 4 1 2]);
    mir.PMI = sum(mir.PMI, [3,4]);
    mir2.PMIraw = permute(mir2.PMIraw, [3 1 2]);
    mir2.PMIraw = sum(mir2.PMIraw, [2,3])
    for i=1:length(mir.algorithms) %number of algorithms
        mir.mir(i,:) = mir.PMI(i,:)./mir2.PMIraw'*100;
    end
end;
tmpname = cell(1,length(allres)); tmpname(:) = { 'name' };
algnames = cellfun(@getfield, { allres.alginfo }, tmpname, 'UniformOutput', false);
for algo=1:length(newalgs)
    inda   = strmatch(lower(newalgs{algo}), lower(mir.algorithms), 'exact');
    indb   = strmatch(lower(newalgs{algo}), lower(pmi.algorithms), 'exact');
    indres = strmatch(lower(newalgs{algo}), lower({ allres.legend }), 'exact');
   
    allmir(:,algo) = mir.mir(inda,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    allres(indres).mir = allmir(:,algo);
    if isfield(mir, 'std')
        allmirstd(:,algo) = mir.std(inda,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    end;
    
    inda = strmatch(newalgs{algo}, algnames, 'exact');
    alldip4(algo) = allres(inda).rvbelow4/length(datasetrange)/71*100;
    alldip5(algo) = allres(inda).rvbelow5/length(datasetrange)/71*100;
    alldip6(algo) = allres(inda).rvbelow6/length(datasetrange)/71*100;
    alldip10(algo) = allres(inda).rvbelow10/length(datasetrange)/71*100;
    alldip15(algo) = allres(inda).rvbelow15/length(datasetrange)/71*100;
    alldip75(algo) = allres(inda).rvbelow75/length(datasetrange)/71*100;
    alldip4subj(algo,:)  = allres(inda).rvbelow4subj /71*100;
    alldip5subj(algo,:)  = allres(inda).rvbelow5subj /71*100;
    alldip6subj(algo,:)  = allres(inda).rvbelow6subj /71*100;
    alldip10subj(algo,:) = allres(inda).rvbelow10subj/71*100;
    alldipp(algo,:) = allres(inda).rvbelowp/length(datasetrange)/71*100';
end;

figure('position', [10 10 900 900]); hold on;
    alpha = 18% strmatch(lower('Picard'), lower(mir.algorithms), 'exact');
    beta = 1 %strmatch(lower('Amica'), lower(mir.algorithms), 'exact');
    gamma = 6 %strmatch(lower('FastICA'), lower(mir.algorithms), 'exact');

    [sortedAmica,sortIdx] = sort(allmir(:,beta),'ascend')
    
    hold on
    alphaArr = allmir(:,alpha);
    
    plot(alphaArr(sortIdx)/1000,'.' ,'color', 'r', 'markersize', 25);
    qw{2} = plot(alphaArr(sortIdx)/1000,'color', 'r', 'LineWidth', 1.5);
    newalgs{alpha}


    betaArr = allmir(:,beta);
    plot(betaArr(sortIdx)/1000, '.' ,'color', [0 0 1], 'markersize', 25);
    qw{1} = plot(betaArr(sortIdx)/1000,'color', [0 0 1], 'LineWidth', 1.5);
    newalgs{beta}



    gammaArr = allmir(:,gamma);
    plot(gammaArr(sortIdx)/1000, '.' ,'color', 'k', 'markersize', 25);
    qw{3} = plot(gammaArr(sortIdx)/1000, 'color', 'k', 'LineWidth', 1.5);
    newalgs{gamma}
    deltaArr = allmir(:,20);
    plot(deltaArr(sortIdx)/1000, '.' ,'color', [0 0.5 0], 'markersize', 25);
    qw{4} = plot(deltaArr(sortIdx)/1000, 'color', [0 0.5 0], 'LineWidth', 1.5);
    xlabel('Dataset');
    ylabel('Mutual information reduction (kbits/sec)');
    setfont(gcf, 'fontsize', 20);
    legend([qw{:}], {'Amica','\color{red} Picard','FastICA', 'PCA'}, 'location', 'best')
    xticklabels = [1 2 3 4 5 6 7 8 9 10 11 12 13];
    xticklabelsSort = xticklabels(sortIdx);
    set(gca,'xtick',[1 2 3 4 5 6 7 8 9 10 11 12 13],'xticklabel',xticklabelsSort)
    yticks('auto')
    yticklabels('auto')

  



% plot all the possible slope and r^2
% -----------------------------------
allr2    = zeros(1,100);
allslope = zeros(1,100);
allalpha = zeros(1,100);
for percent = 1:100
    [ypred allalpha(percent) allr2(percent) allslope(percent)] = myregress(mean(allmir), alldipp(:,percent));
end;



fig_mir_vs_ND = gcf;

figure('position', [2 2 900 900]); 
yyaxis right
plot(-log10(allalpha),'--','LineWidth',2);
xlabel('RV threshold (%)');
ylabel('significance (10^{-X})');
setfont(gcf, 'fontsize', 20);
yyaxis left
%figure('position', [1015 262 277 159]);
plot(allr2,'LineWidth',2);
xlabel('Near dipolar cutoff (% r.v.)');
ylabel('R^{2} value');
setfont(gcf, 'fontsize', 20);
xline(5,'--','LineWidth',2)
legend('R^{2}','','significance')


fig_r2 = gcf;

if plotmir
    saveas(gcf,'MIR_vs_ND_R2.png')
else
    saveas(gcf,'PMI_vs_ND_R2.png')
end


figure('position', [2 2 900 900]); hold on;
amicaIdx = strmatch('Amica', newalgs, 'exact');
for algo=1:length(newalgs)
    amicaMirArr(algo) = abs(mean(allmir(:,algo))/1000 - mean(allmir(:,amicaIdx))/1000);
    mirArr(algo) = mean(allmir(:,algo))/1000;
    %mirArrStd(algo) = std(allmir(:,algo)/1000);
end;
% sort A in descending order (decreasing A values) 
% and keep the sort index in "sortIdx"
[mirArr,sortIdx] = sort(mirArr,'descend');
% sort B using the sorting index
amicaMirArr = amicaMirArr(sortIdx);
yyaxis left
plot(1:length(newalgs), amicaMirArr,'-+','color','k','LineWidth',1.5,'MarkerSize',8)
ylim([min(amicaMirArr) max(amicaMirArr)])
ylabel('MIR difference (\Delta kbits/sec)')
yyaxis right
ax = gca();
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
plot(1:length(newalgs), mirArr,'-o','color','k','LineWidth',1.5,'MarkerSize',8)

picardNewAlgIdx = find(strcmp(newalgs(sortIdx),'Picard'));
picardONewAlgIdx = find(strcmp(newalgs(sortIdx),'Picard-O'));
plot(picardNewAlgIdx, mirArr(picardNewAlgIdx),'.','color','r','MarkerSize',40) %TODO: remove the magic number 4
plot(picardONewAlgIdx, mirArr(picardONewAlgIdx),'.','color','r','MarkerSize',40)

set ( gca, 'ydir', 'reverse' )
ylabel('Mutual information reduction (kbits/sec)')
set(gca,'xtick',[1:length(newalgs)],'xticklabel',newalgsColor(sortIdx))
ylim([min(mirArr) max(mirArr)])
set(gca,'ytick',linspace(min(mirArr),max(mirArr),5)) 
setfont(gcf, 'fontsize', 20);
ytickformat('%,.2f')
saveas(gcf,'MI_difference.png')
print('MI_difference','-depsc')
if 0
% plot MIR or PMI for each subject versus number of dipoles
% --------------------------------
%newalgs = { 'Amica'    'Infomax'    'Ext. Infomax'    'Pearson'    'SHIBBS'    'JADE'    'FastICA' 'TICA'    'JADE opt.'    'SOBI'    'JADE-TD'    'SOBIRO'  };
alldip4subj = alldip4subj';
alldip5subj = alldip5subj';
alldip6subj = alldip6subj';
%%
figure; 
plot2d = 0;
for index = 1:size(allmir,1)
    allvalx = allmir(index,:); %-allmir(index,end);
    allvalx = [ allvalx(end) allvalx(1:end-1) ]; % place sphering before AMICA
    allvalz = allpmi(index,:); %-allmir(index,end);
    allvalz = [ allvalz(end) allvalz(1:end-1) ];
    allvaly = alldip5subj(index,:);
    allvaly = [ allvaly(end) allvaly(1:end-1) ];
    sizemarker = (alldip6subj(index,:)-alldip4subj(index,:))+5;
    %[tmr neworder] = sort(allvalx);
    %allvalx = allvalx(neworder);
    %allvaly = allvaly(neworder);
    indcolor = index;
    if indcolor == 8, indcolor = 13; elseif indcolor ==13, indcolor = 8; end;
    if plot2d == 2
        allvalx = allvalx(2:end) - allvalx(end);
        allvaly = allvaly(2:end) - allvaly(end);
        for pnts = 1:length(allvalx)
            plot( allvalx(pnts), allvaly(pnts), '.', 'markersize', sizemarker(pnts), 'color', allcolors(indcolor,:)); hold on;
        end;
        %plot( allvalx(2:end), allvaly(2:end), '-', 'color', allcolors(indcolor,:));
        [ypred tmp r2(index) slope(index)] = myregress(allvalx, allvaly);
        plot(allvalx, ypred, '-', 'linewidth', 2, 'color', allcolors(indcolor,:));
    elseif plot2d == 1
        allvalx = allvalx - allvalx(end);
        allvaly = allvaly - allvaly(end);
        for pnts = 1:length(allvalx)
            plot( allvalx(pnts), allvaly(pnts), '.', 'markersize', sizemarker(pnts), 'color', allcolors(indcolor,:)); hold on;
        end;
        plot( allvalx(2), allvaly(2), 'x', 'markersize', 10, 'color', allcolors(indcolor,:));
        plot( allvalx(2), allvaly(2), 'o', 'markersize', 10, 'color', allcolors(indcolor,:));
        plot( allvalx(2:end), allvaly(2:end), '-', 'color', allcolors(indcolor,:));
        plot( allvalx(1:2), allvaly(1:2), '--', 'color', allcolors(indcolor,:));

        % plot curve
        if 0
            [ypred] = myregress(allvalx(2:end), allvaly(2:end));
        else
            p = polyfit(allvalx(2:end), allvaly(2:end),2);
            ypred = polyval(p, allvalx(2:end));
        end;
        plot(allvalx(2:end), ypred, '-', 'linewidth', 2, 'color', allcolors(indcolor,:));
    else
        greycol = [0.8 0.8 0.8];
        for pnts = 1:length(allvalx)
            if pnts == 1, marker = 'o'; sizemarker(pnts) = sizemarker(pnts)/3; else marker = '.'; end;
            plot3( allvalx(pnts), allvalz(pnts), allvaly(pnts), marker, 'markersize', sizemarker(pnts), 'color', allcolors(indcolor,:)); hold on;
            plot3( allvalx(pnts),           0.6, allvaly(pnts), marker, 'markersize', sizemarker(pnts), 'color', greycol); hold on;
            %plot3( allvalx(pnts), allvalz(pnts),             0, marker, 'markersize', sizemarker(pnts), 'color', greycol); hold on;
            plot3(        60000, allvalz(pnts), allvaly(pnts), marker, 'markersize', sizemarker(pnts), 'color', greycol); hold on;
        end;
        l = length(allvalx)-1;
        plot3( allvalx(2), allvalz(2), allvaly(2), 'x', 'markersize', 10, 'color', allcolors(indcolor,:)); hold on;
        plot3( allvalx(2), allvalz(2), allvaly(2), 'o', 'markersize', 10, 'color', allcolors(indcolor,:));
        plot3( allvalx(2:end), allvalz(2:end), allvaly(2:end), '-', 'color', allcolors(indcolor,:));
        plot3( allvalx(2:end),  ones(1,l)*0.6, allvaly(2:end), '-', 'color', greycol);
        plot3(ones(1,l)*60000, allvalz(2:end), allvaly(2:end), '-', 'color', greycol);
        plot3( allvalx(1:2), allvalz(1:2), allvaly(1:2), '--', 'color', allcolors(indcolor,:));
        plot3( allvalx(1:2),    [0.6 0.6], allvaly(1:2), '--', 'color', greycol);
        plot3([60000 60000], allvalz(1:2), allvaly(1:2), '--', 'color', greycol);

        % plot curve
%         if 0
%             [ypred] = myregress(allvalx(2:end), allvaly(2:end));
%         else
%             p = polyfit(allvalx(2:end), allvaly(2:end),2);
%             ypred = polyval(p, allvalx(2:end));
%         end;
%         plot(allvalx(2:end), ypred, '-', 'linewidth', 2, 'color', allcolors(indcolor,:));
        xlabel('MIR');
        ylabel('PMI');
        ylim([0.2 0.6])
        zlabel('% of dipoles');
        grid on;
        rotate3d;
        view(-21.5, 36);
        
        meanmir(index) = mean(allvalx);
        meanpmi(index) = mean(allvalz);
    end;
end;

if ~plot2d  
    [ypred alpha r] = myregress(meanmir, meanpmi);
    hold on; plot3(meanmir, ypred, zeros(1,length(ypred)), 'k', 'linewidth', 1);
end;
end


% ----------------------------------------------------------------
% ----------------------------------------------------------------

%return;


%return

if 1
% show numerical results for table
% --------------------------------
%algorithms = {  'Amica' 'Ext. Infomax' 'Pearson' 'Infomax' 'ERICA' 'SONS' ...
%                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
%                'eeA'  'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
%                'AMUSE'  'Sphering' 'PCA'    }; % 'icaML' 'promax' 'unica'    'simbec' 
algorithms = newalgs
for algo=1:length(algorithms)
    inda = strmatch(algorithms{algo}, mir.algorithms, 'exact');
    allmir(:,algo) = (mir.mir(end,datasetrange)-mir.mir(inda,datasetrange))*1.4427;
    inda = strmatch(algorithms{algo}, algnames, 'exact');
    if ~isempty(inda)
         alldip10table(algo) = allres(inda).rvbelow10/length(datasetrange)/71*100;
    else alldip10table(algo) = 0;
    end;
end;
[tmp sorti] = sort(mean(allmir,1));
allmir    = allmir(:,sorti);
newalgs2  = algorithms(sorti);
alldip10table = alldip10table(sorti);
[tmp sorti] = sort(mean(allmir(:,1:end-2),2));
allmir    = allmir(sorti,:);

meanmir = mattocell(mean(allmir,1));
stdmir  = mattocell(std(allmir,[],1));
alldip10tablecell = mattocell(alldip10table);
%{ newalgs2{:}; meanmir{:}; alldip10tablecell{:} }'
for i = 1:length(meanmir)
    fprintf('%10s\t%2.2f\t%2.1f\n', newalgs2{i}, meanmir{i}, alldip10tablecell{i})
end;

%return
end

