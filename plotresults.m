% This script plots results for all ICA algorithm
% Because some ICA algorithm are computed and not plotted
% you must add the name of your algorithm to the list
% Arnaud Delorme - Aug 2011

clear;
newalgs = {  'AMICA' 'Ext. Infomax' 'Pearson' 'Infomax' ...
                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
                 'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
                'AMUSE'  'PCA'  'Sphering' }; % 'icaML' 'promax' 'unica'    'simbec' 

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
        allres(countplot).legend   = allalgs(ALGONUM).name;
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
    allres(index).rvbelow4    = length(find(allres(index).rvalgo(:) < 0.045));
    allres(index).rvbelow5    = length(find(allres(index).rvalgo(:) < 0.05));
    allres(index).rvbelow6    = length(find(allres(index).rvalgo(:) < 0.055));
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
figure('position', [2   262   728   662]);
hh = semilogy(allres(1).allrvsorted, 'w'); hold on;
delete(hh);
for index = 1:length(allres)
    %tmpind = find( allres.rvalgo(:,index) );
    %tmpdat = allres.rvalgo(tmpind,index);
    %tmpind = ~isnan(tmpdat);
    %tmpdat = tmpdat(tmpind);
    %abscicia = linspace(1,72, length(tmpdat));
    %if length(tmpdat) ~= 0
    %   h = semilogy(abscicia, tmpdat); hold on;
    h = semilogy(linspace(1,100,length(allres(index).allrvsorted)), allres(index).allrvsorted); hold on;
    %h = plot([1:inc:71.8], log10(allres.rvalgo(:,index))); hold on;
    set(h, 'linewidth', 2, 'color', allcolors(mod(index-1,13)+1,:));
    if index>12, set(h, 'linestyle', '--'); end;
end;
%figure; plot(sqrt([1:71*5]), log10(allres.rvalgo'));
%figure; semilogx(log10(allres.rvalgo'));
hh = semilogy(linspace(1,100, 710), nan_mean(nan_mean(totallrv,2),2), 'k-.', 'linewidth', 4); hold on;
lh = legend({ allres.legend 'raw EEG' });
view([90 270]);
xlabel('Percent of ICA components');
ylabel('Dipole model residual variance (%)');
set(gca, 'ylim', [0.004, 1]);
set(gca, 'xlim', [0 100]);
set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
set(gcf, 'paperpositionmode', 'auto');
setfont(gcf, 'fontsize', 16);
set(lh, 'position', [0.1598 0.2119 0.1670 0.7141]);
%eval(sprintf('print -depsc icacomp%d-%d.eps', datasetrange(1), datasetrange(end)));
%print -djpeg icacomp.jpg

% plot PMI versus MIR
% -------------------
clear allmir allpmi;
mir = load( '-mat', 'pmi_save.mat');
mir2 = load('-mat', 'pmi_save_data.mat');
mir.mir = mir.PMI./repmat(mir2.PMI, [size(mir.PMI,1) 1])*100;
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
pmi = load( '-mat', 'pmi_stat_save_oct2011.mat');
pmi2 = load('-mat', 'pmi_stat_save_raw_oct2011.mat');
%pmi.pmi = pmi.PMI*1.4427*250*2845/1000;
pmi.pmi = pmi.PMI./repmat(pmi2.PMIraw, [size(pmi.PMI,1) 1])*100; % in %
if plotmir % CHANGE HERE FOR PMI VERSUS MIR
    mir = load('-mat', 'mir_new.mat');
    mir.mir = mir.mir*1.4427*250;
else
    %mir = load( '-mat', 'pmi_save.mat');
    %mir2 = load('-mat', 'pmi_save_data.mat');
    %mir.mir = mir.PMI./repmat(mir2.PMI, [size(mir.PMI,1) 1])*100;
    mir = load( '-mat', 'pmi_stat_save_oct2011.mat');
    mir2 = load('-mat', 'pmi_stat_save_raw_oct2011.mat');
    mir.mir = mir.PMI./repmat(mir2.PMIraw, [size(mir.PMI,1) 1])*100;
    mir.std = 3*sqrt(mir.vPMI./repmat(mir2.PMIraw, [size(mir.vPMI,1) 1]))*100;
end;
tmpname = cell(1,length(allres)); tmpname(:) = { 'name' };
algnames = cellfun(@getfield, { allres.alginfo }, tmpname, 'UniformOutput', false);
for algo=1:length(newalgs)
    inda   = strmatch(lower(newalgs{algo}), lower(mir.algorithms), 'exact');
    indb   = strmatch(lower(newalgs{algo}), lower(pmi.algorithms), 'exact');
    indres = strmatch(lower(newalgs{algo}), lower({ allres.legend }), 'exact');
    
    %allmir(:,algo) = (mir.mir(end,datasetrange)-mir.mir(inda,datasetrange))*1.4427; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits | 250 is the sampling rate
    allmir(:,algo) = mir.mir(inda,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    allpmi(:,algo) = pmi.pmi(indb,datasetrange); % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    allpmi(:,algo) = pmi2.PMIraw(datasetrange)*1.4427*250*2845/1000; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    allres(indres).mir = allmir(:,algo);
    allres(indres).pmi = allpmi(:,algo);
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
figure('position', [732 504 560 420]); hold on;
for i=1:length(alldip10), 
    plot(mean(allmir(:,i)), alldip5(i), '.', 'color', 'k', 'markersize', 18); 
    if isfield(mir, 'std')
        valx1 = (mean(allmir(:,i))-mean(allmirstd(:,i)));
        valx2 = (mean(allmir(:,i))+mean(allmirstd(:,i)));
%         plot([valx1 valx1], [ alldip4(i) alldip6(i) ], 'k');
%         plot([valx2 valx2], [ alldip4(i) alldip6(i) ], 'k');
%         plot([valx1 valx2], [ alldip4(i) alldip4(i) ], 'k');
%         plot([valx1 valx2], [ alldip6(i) alldip6(i) ], 'k');
        ellipse(valx2-valx1,alldip6(i)-alldip4(i),0,(valx1+valx2)/2,(alldip4(i)+alldip6(i))/2,'k')
    end;
    %if i < 13
    %     plot(mean(allmir(:,i)), alldip5(i), 'x', 'color', allcolors(i,:)   , 'markersize', 5); 
    %else plot(mean(allmir(:,i)), alldip5(i), 'o', 'color', allcolors(i-12,:), 'markersize', 5); 
    %end;
    if plotmir
        text(mean(allmir(:,i))/1000+0.02, alldip5(i), newalgs{i});
    else
        text(mean(allmir(:,i))/1000+0.00002, alldip5(i), newalgs{i});
    end;
end;

%% xlim([-0.05 2]*1.4427)
[ypred alpha5 r5 slope] = myregress(mean(allmir)/1000, alldip5);
fprintf('-----------------------------------\nRegression at 5%% threshold: r^2=%1.3f\tp=%e\n', r5, alpha5);
hold on; plot(mean(allmir)/1000, ypred, 'k--', 'linewidth', 1.5);
if 0
    %for i=1:length(alldip10), plot(mean(allmir(:,i)), alldip75(i), '.', 'color', allcolors(i,:), 'markersize', 18); end;
    %[ypred alpha75 r75 slope] = myregress(mean(allmir), alldip75);
    %hold on; plot(mean(allmir), ypred, 'k', 'linewidth', 1.5);
    
    for i=1:length(alldip10), plot(mean(allmir(:,i)), alldip10(i), '.', 'color', allcolors(i,:), 'markersize', 18); end;
    [ypred alpha10 r10 slope] = myregress(mean(allmir), alldip10);
    hold on; plot(mean(allmir), ypred, 'k:', 'linewidth', 2);
end;
setfont(gcf, 'fontsize', 16);
if plotmir
    xlabel('Mutual information reduction (kbits/sec)');
else
    xlabel('Remaining Pairwise Mutual Information  (10^-3 %)', 'interpreter', 'none');
end;
ylabel('Percent near-dipolar (r.v. < 5%) components');

% plot all the possible slope and r^2
% -----------------------------------
allr2    = zeros(1,100);
allslope = zeros(1,100);
allalpha = zeros(1,100);
for percent = 1:100
    [ypred allalpha(percent) allr2(percent) allslope(percent)] = myregress(mean(allmir), alldipp(:,percent));
end;

figure('position', [732 262 281 159]); plot(-log10(allalpha));
xlabel('RV threshold (%)');
ylabel('significance (10^{-X})');
setfont(gcf, 'fontsize', 16);

figure('position', [1015 262 277 159]); plot(allr2);
xlabel('RV threshold (%)');
ylabel('R^{2} value');
setfont(gcf, 'fontsize', 16);

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

%% plot dataset number 8 only
% --------------------------
figure; 
index = 8;
allvalx = allmir(index,:); %-allmir(index,end);
allvaly = alldip5subj(index,:);
[tmr neworder] = sort(allvalx);
allvalxx = allvalx(neworder);
allvalyy = allvaly(neworder);
plot( allvalxx, allvalyy, '.', 'markersize', 10, 'color', allcolors(index,:));
%plot( allvalx, allvaly, '-', 'color', allcolors(index,:));
hold on;

% plot curve
[ypred] = myregress(allvalxx, allvalyy);
plot(allvalxx, ypred, '-', 'linewidth', 2, 'color', allcolors(index,:));
for index = 1:size(allmir,2)
    text( allvalx(index)+30, allvaly(index), newalgs{index});
end;

% ----------------------------------------------------------------
% ----------------------------------------------------------------

return;

figure; 
algo = 2
plot(allmir(:,algo), alldip5subj(algo,:), '.');
[ypred alpha r slope] = myregress(allmir(:,algo), alldip5subj(algo,:));
hold on; plot(allmir(:,algo), ypred);
title(sprintf('AMICA: Alpha=%2.0f%% (r2=%0.2f)', alpha*100, r*100)); 
xlabel('Mutual information reduction (kbits/sec)');
ylabel('Number of near-dipolar comps. (%)');
setfont(gcf, 'fontsize', 16)

figure; plot(allslope*1000);
xlabel('RV threshold (%)');
ylabel('Slope value');
setfont(gcf, 'fontsize', 16);

% plot single dataset number of dipoles for AMICA
% -----------------------------------------------
for index = 1:length(allres)
    for s = 1:size(allres(index).rvalgo,2)
        allres(index).rvbelow5subj(s) = length(find(allres(index).rvalgo(:,s) < 0.05));
        allres(index).rvbelow10subj(s) = length(find(allres(index).rvalgo(:,s) < 0.1));
    end;
end;
mir = load('-mat', 'mir_new.mat');
for algo=1:length(newalgs)
    inda = strmatch(newalgs{algo}, mir.algorithms, 'exact');
    allmir(:,algo) = (mir.mir(end,datasetrange)-mir.mir(inda,datasetrange))*1.4427; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits | 250 is the sampling rate
    allmir(:,algo) = mir.mir(inda,datasetrange)*1.4427*250; % 1 Nats = 2*log2(exp(1)) = 1.4427 bits
    alldip5subj(algo,:)  = allres(algo).rvbelow5subj /71*100;
    alldip10subj(algo,:) = allres(algo).rvbelow10subj/71*100;
end;
figure; 
algo = 2
plot(allmir(:,algo), alldip5subj(algo,:), '.');
[ypred alpha r slope] = myregress(allmir(:,algo), alldip5subj(algo,:));
hold on; plot(allmir(:,algo), ypred);
title(sprintf('AMICA: Alpha=%2.0f%% (r2=%0.2f)', alpha*100, r*100)); 
xlabel('Mutual information reduction (kbits/sec)');
ylabel('Number of near-dipolar comps. (%)');
setfont(gcf, 'fontsize', 16)

return

% show numerical results for table
% --------------------------------
algorithms = {  'Amica' 'Ext. Infomax' 'Pearson' 'Infomax' 'ERICA' 'SONS' ...
                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
                'eeA'  'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
                'AMUSE'  'Sphering' 'PCA'    }; % 'icaML' 'promax' 'unica'    'simbec' 
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
return

% compute difference with EEG data
% --------------------------------
allrveeg   = reshape(nan_mean(tmpallrv,2), [71 size(tmpallrv,3)]);
percentout = 10.^[0:0.025:2];
clear allrveeg2;
for index = 1:length(datasetrange)
    nonnanindices(index,:) = setdiff([1:71], find(isnan(allrveeg(:,index))));
    allrveeg2(index,:)     = vectdata(nonnanindices(index,:), ...
                             100*allrveeg(nonnanindices(index,:),index)', 'timesout', percentout);
end;

%figure; % testing interpolation
%subplot(1,2,1); semilogx(percentout, allrveeg2');
%subplot(1,2,2); semilogy([1:71], allrveeg); view([90 270]);

% plot result of the algorithms against results of the EEG
% interpolating for residual variance 1:100 in log space
% ------------------------------------------------------
figure;
allres.allrvalgo = sort(allres.allrvalgo,2);
for index = 1:size(allres.allrvalgo,1) % scan algos
    
    % compute difference for algo
    warning off MATLAB:griddata:DuplicateDataPoints
    clear tmpres;
    for dataset = 1:length(datasetrange)
        tmpnonnan = setdiff([1:71], find(isnan(allres.allrvalgo(index,:,dataset))));
        tmpres(dataset,:) = vectdata(tmpnonnan, ...
                             100*allres.allrvalgo(index,tmpnonnan,dataset), 'timesout', percentout);
    end;
    tmpres = nan_mean(tmpres - allrveeg2,1);
    h = semilogx( percentout, tmpres); hold on;
    if index > 13
        set(h, 'linewidth', 2, 'linestyle', '--', 'color', allcolors(index-13,:));
    else
        set(h, 'linewidth', 2, 'linestyle', '-',  'color', allcolors(index,:));
    end;
end;
plot([0.01 100], [0 0], 'k', 'linewidth', 2);
lh = legend(allres.legend{:});
xlabel('Residual variance (%)');
ylabel('Number of components compared to EEG');
set(gca, 'xlim', [0.4, 100]);
set(gca, 'ylim', [-60 20]);
set(gca, 'xtick', [1 5 10 20 100], 'xticklabel', [1 5 10 20 100]);
set(gcf, 'paperpositionmode', 'auto');
set(gcf, 'position', [245.0000  302.0000  843.0000  793.0000]);
title(sprintf('Datasets %d-%d', datasetrange(1), datasetrange(end)));
setfont(gcf, 'fontsize', 16);
set(lh, 'position', [0.1444    0.2942    0.1756    0.6142]);
eval(sprintf('print -depsc icacompdiff%d-%d.eps', datasetrange(1), datasetrange(end)));
%print -djpeg icacomp.jpg

% plot results for runica 10 datasets
% ------------------------------------------------------
allres.allrvalgo = sort(allres.allrvalgo,2);
for index = 1:6 % scan algos
    figure;
    
    % compute difference for algo
    warning off MATLAB:griddata:DuplicateDataPoints
    clear tmpres;
    for dataset = 1:length(datasetrange)
        tmpnonnan = setdiff([1:71], find(isnan(allres.allrvalgo(index,:,dataset))));
        tmpres    = vectdata(tmpnonnan, ...
                             100*allres.allrvalgo(index,tmpnonnan,dataset), 'timesout', percentout);
        tmpres = nan_mean(tmpres - allrveeg2(dataset,:),1);
        h = semilogx( percentout, tmpres); hold on;
        alllh{dataset} = [ 'Dataset ' int2str(dataset) ];
        if dataset > 13
            set(h, 'linewidth', 2, 'linestyle', '--', 'color', allcolors(dataset-13,:));
        else
            set(h, 'linewidth', 2, 'linestyle', '-',  'color', allcolors(dataset,:));
        end;
    end;
    plot([0.01 100], [0 0], 'k', 'linewidth', 2);
    lh = legend(alllh);
    xlabel('Residual variance (%)');
    ylabel('Number of components compared to EEG');
    set(gca, 'xlim', [0.4, 100]);
    set(gca, 'ylim', [-35 35]);
    set(gca, 'xtick', [1 5 10 20 100], 'xticklabel', [1 5 10 20 100]);
    set(gcf, 'paperpositionmode', 'auto');
    set(gcf, 'position', [245.0000  302.0000  843.0000  793.0000]);
    title(sprintf('Datasets %d-%d', datasetrange(1), datasetrange(end)));
    setfont(gcf, 'fontsize', 16);
    set(lh, 'position', [0.1444    0.2942    0.1756    0.6142]);
    title(allres.legend{index});
    eval(sprintf('print -depsc icadiff%d.eps', index));
end;
%print -djpeg icacomp.jpg

% plot algorithm against time
% ---------------------------
figure;
hh = semilogx(allres.algtime, allres.rvbelow10/355, '.'); hold on;
delete(hh);
for index = 1:size(allres.rvalgo,2) 
    %hh = errorbar(allres.algtime(index), rvbelow10p(index)/5, allres.stdrvbelow10(index), allres.stdrvbelow10(index)); hold on;
    hh        = plot(allres.algtime(index), allres(index).rvbelow10p(index), '.'); hold on;
    hh(end+1) = plot([allres(index).algtime allres(index).algtime], [allres(index).rvbelow10p-allres(index).stdrvbelow10 allres(index).rvbelow10p+allres(index).stdrvbelow10]);
    hh(end+1) = plot([allres(index).algtime*1.05 allres(index).algtime/1.05], [allres(index).rvbelow10p(index)-allres(index).stdrvbelow10 allres(index).rvbelow10p-allres(index).stdrvbelow10]);
    hh(end+1) = plot([allres(index).algtime*1.05 allres(index).algtime/1.05], [allres(index).rvbelow10p(index)+allres(index).stdrvbelow10 allres(index).rvbelow10p+allres(index).stdrvbelow10]);
    set(hh, 'color', allcolors(mod(index-1,12)+1,:), 'markersize', 12);
    text( allres(index).algtime*0.85, allres(index).rvbelow10p+1, allres(index).legend);
end;
view([90 270]);
title(sprintf('Datasets %d-%d', datasetrange(1), datasetrange(end)));
set(gca, 'ylim', [0 50], 'xticklabel', [10 100 1000 10000 100000 1000000]);
ylabel('Number of dipolar components (rv < 10%)');
xlabel('Approximate computation time (s)');
set(gcf, 'paperpositionmode', 'auto');
eval(sprintf('print -depsc icatime%d-%d.eps', datasetrange(1), datasetrange(end)));
%print -depsc icatime.eps
%print -djpeg icatime.jpg

% correlation of the best algorithms
% ----------------------------------
allwinv = reshape(winv, [24 71 71 14]);
allwinv = permute(allwinv([1:10,12:14],:,:,:), [2 3 1 4]);
clear allcorrs allwinv2;
for DATASET = 1:14
    [allcorrs(:,:,DATASET) allwinv2(:,:,:, DATASET)] = matcorrm(allwinv(:,:,:, DATASET), allwinv(:,:,1, DATASET), ...
                                                      'orient', 'rows', 'recur', 'on');
end;

plotdensities;

% plot correlated images
% ----------------------
figure;
for index = 1:13
    subplot(3,5,index); topoplot(allwinv2(:,1,index,1), EEG.chanlocs, 'electrodes', 'off');
end;

mean(allcorrs(:))
corr1  = mean(mean(allcorrs(1 ,:,:),2),3);
corr35 = mean(mean(allcorrs(35,:,:),2),3);
corr71 = mean(mean(allcorrs(71,:,:),2),3);
    

for index = 1:32
    maxc(index) = max(EEG.icawinv(index,:));
end;
sort(maxc)
