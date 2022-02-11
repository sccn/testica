% find dipoles in EEG data, reconstructed data from ICA/PCA (shuffled activity)
% have to define 1 values
% DATASET: sample dataset to use
% ALGO: 1=ica  or 2=pca
% TYPE: 1=data or 2=ica/pca

compute = 0;
plot    = 1;

if compute == 1
    if exist('DATASET') == 1
        switch DATASET
         case 1 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/km81/Sternberg/');
         case 2 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/jo74/Sternberg/');
         case 3 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ds76/Sternberg/');
         case 4 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/cj82/Sternberg/');
         case 5 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ap82/Sternberg/');
         case 6 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ke70/Sternberg/');
         case 7 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/tp62/Sternberg/');
         case 8 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/cz84/Sternberg/');
         case 9 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/gm84/Sternberg/');
         case 10, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/gv84/Sternberg/');
         case 11, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/nf68/Sternberg/');
         case 12, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ds80/Sternberg/');
         case 13, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/kb77/Sternberg/');
         case 14, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ts79/Sternberg/');
        end;
        EEG = eeg_checkset(EEG);
        EEG = pop_select(EEG, 'trial', setdiff([1:EEG.trials], [1:4:EEG.trials]));
    end;
    load('-mat', [ '../ica' int2str(DATASET) '_72_02_runica.mat' ]); % load runica by default (for subtracting components
    EEG.icaweights = W;
    EEG.icasphere  = eye(size(W,2));
    EEG.icawinv    = inv(W);
    EEG.icaact     = [];
    EEG = eeg_checkset(EEG);
    
    if ALGO == 1
        FILENAME = [ '../ica' int2str(DATASET) '_72_02_runica.mat' ];
        ALGONAME = 'runica';
    elseif ALGO == 2
        FILENAME = [ '../ica' int2str(DATASET) '_72_43_pca.mat' ];
        ALGONAME = 'pca';
    else
        FILENAME = [ '../ica' int2str(DATASET) '_72_08_jader.mat' ];
        ALGONAME = 'jader';
    end;
    
    tmppwd = pwd;
    cd ..;
    compinfo;
    cd(tmppwd);
    if TYPE == 1
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data';
    elseif TYPE == 2
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data_shuffled';
    elseif TYPE == 3
        EEG = pop_subcomp(EEG, eyes{DATASET});
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data_noeye';
    elseif TYPE == 4
        EEG = pop_subcomp(EEG, [ eyes{DATASET} muscle{DATASET} ]);
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data_noeyemuscle';
    elseif TYPE == 5
        EEG = pop_subcomp(EEG, [ eyes{DATASET} muscle{DATASET} singchan{DATASET} dccomp{DATASET}]);
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data_noartifacts';
    elseif TYPE == 6
        EEG = pop_subcomp(EEG, setdiff([1:71], bio{DATASET}));
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.data, EEG.chanlocs, 'nfits', 10);
        ALGONAME = 'data_bio';
    else        
        load('-mat', FILENAME);
        EEG.icaweights = W;
        EEG.icasphere  = eye(size(W,2));
        EEG.icawinv    = inv(W);
        EEG = eeg_checkset(EEG);
        [cumout,lowrvs,posxyz,momxyz,idx] = spatcomp(EEG.icaact,EEG.chanlocs, 'nfits', 10, ...
                                                     'icamaps', EEG.icawinv, 'icamix', 'on');
    end;
    filename = sprintf('randomdip_%d_%s.mat', DATASET, ALGONAME);
    save(filename, '-mat','cumout', 'lowrvs', 'posxyz', 'momxyz', 'idx');
end;
    
% plot results
% ------------
if plot == 1

    if ~exist('TYPE'), TYPE = 1; end;
    if TYPE == 1,     ALGONAME = 'data';
    elseif TYPE == 2, ALGONAME = 'data_shuffled';
    elseif TYPE == 3, ALGONAME = 'data_noeye';
    elseif TYPE == 4, ALGONAME = 'data_noeyemuscle';
    elseif TYPE == 5, ALGONAME = 'data_noartifacts';
    elseif TYPE == 6, ALGONAME = 'data_bio';
    elseif ALGO == 1, ALGONAME = 'runica';
    elseif ALGO == 2, ALGONAME = 'pca';
    else              ALGONAME = 'jader';
    end;
    
    % add EEG in the background
    % -------------------------
    clear totallrv
    count = 1;
    for index = 1:14
        filename2 = sprintf('randdip/randomdip_%d_%s.mat',index, ALGONAME);
        load('-mat', filename2);
        totallrv(count,:) = cumout(:)';
        count = count+1;
    end;
    totallrv = totallrv';
    tmpallrv = reshape(totallrv, 71, 10, size(totallrv,2));
    totallrv = sort(totallrv,1);
    
    % plot only EEG (10 datasets)
    % ---------------------------
    figure;
    allcolors = hsv(12);
    allcolors(1 ,:) = [1.0000    0.6855    0.1109];
    allcolors(2 ,:) = [0.5 0.5 0.5];
    allcolors(4 ,:) = [0   0.6 0];
    allcolors(6 ,:) = [0.8 0.8 0];
    allcolors(9 ,:) = [0   0   0];
    allcolors(13,:) = [0.7 0   0.7];
    allcolors(14, :) = [0.7 1 0.7];
    hh = semilogy(rand(1,100), 'w'); hold on;
    X  = [1:71];
    for index = 1:size(tmpallrv,3)    
        allmeans(:, index) = nan_mean(tmpallrv(:,:,index),2);
        topy = nan_mean(tmpallrv(:,:,index),2)+nan_std(tmpallrv(:,:,index)')';
        boty = nan_mean(tmpallrv(:,:,index),2)-nan_std(tmpallrv(:,:,index)')';
        hold on; h = fillcurves(X, topy', boty');
        set(h, 'edgecolor', 'none', 'facecolor',  allcolors(index, :), ...
               'facealpha', 0.4);
    end;
    delete(hh);
    view([90 270]);
    xlabel('Number of components');
    ylabel('Residual variance (%)');
    set(gca, 'ylim', [0.01, 1]);
    set(gca, 'xlim', [0 71]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    set(gcf, 'paperpositionmode', 'auto');
    title(ALGONAME, 'interpreter', 'none');
    
    axes('position', [0 0.6 0.4 0.4]);    
    hh = semilogy(rand(1,100), 'w'); hold on;
    errorbar(nan_mean(allmeans,2), nan_std(allmeans')'); hold on;
    view([90 270]);
    set(gca, 'ylim', [0.01, 1]);
    set(gca, 'xlim', [0 71]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    delete(hh);
    
    eval([ 'print -deps rand_' ALGONAME '.eps' ]);
    return;
    
    cumout1 = [];
    cumout2 = [];
    cumout3 = [];
    for DATASET = 1:14
        load('-mat', ['randomdip_' int2str(DATASET) '_runica.mat']);
        cumout1 = [cumout1 cumout];
        load('-mat', ['randomdip_' int2str(DATASET) '_data.mat']);
        cumout2 = [cumout2 cumout];
        load('-mat', ['randomdip_' int2str(DATASET) '_pca.mat']);
        cumout3 = [cumout3 cumout];
    end;
    std1 = nan_std(cumout1')'; cumout1 = nan_mean(cumout1,2);
    std2 = nan_std(cumout2')'; cumout2 = nan_mean(cumout2,2);
    std3 = nan_std(cumout3')'; cumout3 = nan_mean(cumout3,2);
    figure;
    view([90 270]);
    h1 = semilogy(cumout1, 'b'); hold on;
    h2 = semilogy(cumout2, 'r');
    h3 = semilogy(cumout3, 'g');
    view([90 270]);
    errorbar(cumout1, std1/2);
    errorbar(cumout2, std2/2, 'r');
    errorbar(cumout3, std3/2, 'g');
    view([90 270]);
    set(gca, 'ylim', [0.04, 1]);
    set(gca, 'xlim', [0 71]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    xlabel('Number of components');
    ylabel('Residual variance (%)');
    print -depsc randomdip_all2.eps
    
    % plot mean for each subject
    % --------------------------
    cumout1 = [];
    cumout2 = [];
    cumout3 = [];
    for DATASET = 1:14
        load('-mat', ['randomdip_' int2str(DATASET) '_runica.mat']);
        cumout1 = [cumout1 nan_mean(cumout,2)];
        load('-mat', ['randomdip_' int2str(DATASET) '_data.mat']);
        cumout2 = [cumout2 nan_mean(cumout,2)];
        load('-mat', ['randomdip_' int2str(DATASET) '_pca.mat']);
        cumout3 = [cumout3 nan_mean(cumout,2)];
    end;
    figure;
    allleg = {};
    for index = 1:14
        h(1:5)   = semilogy(cumout1(:,index), 'b'); hold on;
        h(6:10)  = semilogy(cumout2(:,index), 'b');
        h(11:15) = semilogy(cumout3(:,index), 'b--');
        allleg = { allleg{:} [ 'runica ' int2str(index) ] [ 'data ' int2str(index) ] [ 'pca ' int2str(index) ]};
        set(h(6:10), 'linewidth', 2);
        switch index
         case 1, set(h, 'color', 'b');
         case 2, set(h, 'color', 'r');
         case 3, set(h, 'color', 'g');
         case 4, set(h, 'color', 'k');
         case 5, set(h, 'color', 'm');
        end;
    end;
    legend(allleg{:});
    view([90 270]);
    set(gca, 'ylim', [0.04, 1]); 
    set(gca, 'xlim', [0 71]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    xlabel('Number of components');
    ylabel('Residual variance (%)');
    print -depsc randomdip_allmean.eps
    
    % plot 710 dipoles
    % ----------------
    cumout1 = [];
    cumout2 = [];
    cumout3 = [];
    for DATASET = 1:14
        load('-mat', ['randomdip_' int2str(DATASET) '_runica.mat']);
        cumout1 = [cumout1 sort(cumout(:))];
        load('-mat', ['randomdip_' int2str(DATASET) '_data.mat']);
        cumout2 = [cumout2 sort(cumout(:))];
        load('-mat', ['randomdip_' int2str(DATASET) '_pca.mat']);
        cumout3 = [cumout3 sort(cumout(:))];
    end;
    figure;
    allleg = {};
    for index = 1:14
        h(1:5)   = semilogy(cumout1(:,index), 'b'); hold on;
        h(6:10)  = semilogy(cumout2(:,index), 'b');
        h(11:15) = semilogy(cumout3(:,index), 'b--');
        allleg = { allleg{:} [ 'runica ' int2str(index) ] [ 'data ' int2str(index) ] [ 'pca ' int2str(index) ]};
        set(h(6:10), 'linewidth', 2);
        allc = hsv(14);
        set(h, 'color', allc(index,:));
    end;
    legend(allleg{:});
    view([90 270]);
    set(gca, 'ylim', [0.04, 1]);
    set(gca, 'xlim', [0 710]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    xlabel('Number of components');
    ylabel('Residual variance (%)');
    print -depsc randomdip_710all.eps
    
    std1 = nan_std(cumout1')'; cumout1mean = nan_mean(cumout1,2);
    std2 = nan_std(cumout2')'; cumout2mean = nan_mean(cumout2,2);
    std3 = nan_std(cumout3')'; cumout3mean = nan_mean(cumout3,2);
    figure;
    h1 = semilogy(cumout1mean, 'b'); hold on;
    view([90 270]);
    errorbar(cumout2mean, std2/10, 'r');
    errorbar(cumout1mean, std1/10, 'b');
    errorbar(cumout3mean, std3/10, 'g');
    view([90 270]);
    set(gca, 'ylim', [0.04, 1]);
    set(gca, 'xlim', [0 710]);
    set(gca, 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    xlabel('Number of components');
    ylabel('Residual variance (%)');
    print -depsc randomdip_710mean2.eps
    
    % plot big figure (6 subplots)
    % ----------------------------
    clear plot
    figure;
    subplot(3,3,1); semilogy(cumout1); hold on; semilogy(nan_mean(cumout2,2), 'k', 'linewidth',2);
    set(gca, 'ylim', [0.04, 1], 'xlim', [0 710],'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    view([90 270]); xlabel('Number of components'); ylabel('Residual variance (%)');
    
    subplot(3,3,2); semilogy(cumout2); hold on; semilogy(nan_mean(cumout2,2), 'k', 'linewidth',2);
    set(gca, 'ylim', [0.04, 1], 'xlim', [0 710],'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    view([90 270]); xlabel('Number of components'); ylabel('Residual variance (%)');
    
    subplot(3,3,3); semilogy(cumout3); hold on; semilogy(nan_mean(cumout2,2), 'k', 'linewidth',2);
    set(gca, 'ylim', [0.04, 1], 'xlim', [0 710],'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    view([90 270]); xlabel('Number of components'); ylabel('Residual variance (%)');
    
    subplot(3,3,4); loglog(cumout2, cumout1);
    hold on; plot([0.01 1], [0.01 1], 'k', 'linewidth', 2);

    set(gca, 'ylim', [0.04, 1], 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    xlabel('Residual variance EEG (%)'); ylabel('Residual variance ICA (%)');
    
    subplot(3,3,5); loglog(cumout2, cumout3);
    hold on; plot([0.01 1], [0.01 1], 'k', 'linewidth', 2);
    set(gca, 'ylim', [0.04, 1], 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    xlabel('Residual variance EEG (%)'); ylabel('Residual variance PCA (%)');

    subplot(3,3,6); loglog(cumout1, cumout3);
    hold on; plot([0.01 1], [0.01 1], 'k', 'linewidth', 2);
    set(gca, 'ylim', [0.04, 1], 'ytick', [0.01 0.05 0.1 0.2 1], 'yticklabel', [1 5 10 20 100]);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    xlabel('Residual variance ICA (%)'); ylabel('Residual variance PCA (%)');

    subplot(3,3,7); edges = 10.^(linspace(log10(0.01), 0, 50)); N = histc(cumout1(:), edges);
    h=bar(edges, N, 'histc'); set(gca, 'xscale', 'log'); %hist(cumout1(:), 50);
    tmpobj = get(gca, 'children'); if tmpobj(1) == h, h=tmpobj(2); else h=tmpobj(1); end; delete(h);
    hold on; plot([0.05 1], [0 700], 'k', 'linewidth', 2);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    %set(gca, 'xlim', [0 1], 'xtick', [0:0.2:1], 'xticklabel', [0:20:100]);
    xlabel('ICA residual variance (%)'); ylabel('# of components'); ylim([0 600]);

    subplot(3,3,8); edges = 10.^(linspace(log10(0.01), 0, 50)); N = histc(cumout2(:), edges);
    h=bar(edges, N, 'histc'); set(gca, 'xscale', 'log'); %hist(cumout2(:), 50);
    tmpobj = get(gca, 'children'); if tmpobj(1) == h, h=tmpobj(2); else h=tmpobj(1); end; delete(h);
    hold on; plot([0.05 1], [0 700], 'k', 'linewidth', 2);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    %set(gca, 'xlim', [0 1], 'xtick', [0:0.2:1], 'xticklabel', [0:20:100]);
    xlabel('EEG residual variance (%)'); ylabel('# of components'); ylim([0 600]);

    subplot(3,3,9); edges = 10.^(linspace(log10(0.01), 0, 50)); N = histc(cumout3(:), edges);
    h=bar(edges, N, 'histc'); set(gca, 'xscale', 'log'); %hist(cumout3(:), 50);
    tmpobj = get(gca, 'children'); if tmpobj(1) == h, h=tmpobj(2); else h=tmpobj(1); end; delete(h);
    hold on; plot([0.05 1], [0 700], 'k', 'linewidth', 2);
    set(gca, 'xlim', [0.04, 1], 'xtick', [0.01 0.05 0.1 0.2 1], 'xticklabel', [1 5 10 20 100]);
    %set(gca, 'xlim', [0 1], 'xtick', [0:0.2:1], 'xticklabel', [0:20:100]);
    xlabel('PCA residual variance (%)'); ylabel('# of components'); ylim([0 600]);
    
    set(gcf, 'position', [ 572   368   853   759]);
    set(gcf, 'paperpositionmode', 'auto');
    print -depsc randsummary.eps
end;

return;

% launch on cluser (copy and paste to command line)
% -------------------------------------------------
for DATASET = 1:14
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 1; ALGO = 1; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_data -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 3; ALGO = 1; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET)  '_ica -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 3; ALGO = 2; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET)  '_pca -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
end;

for DATASET = 1:14
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 2; ALGO = 1; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_2 -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 3; ALGO = 1; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_3 -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 4; ALGO = 3; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_4 -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 5; ALGO = 3; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_5 -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
    fprintf(['echo "cd matlab/mica/randdip/; matlab -nodesktop -nosplash -r ''' ...
             'TYPE = 6; ALGO = 3; DATASET = ' int2str(DATASET) ';' ...
             ' randomica; quit;''" | qsub -N rand_' int2str(DATASET) '_6 -o /home/arno/matlab/mica/randdip/qsub -p -4\n'], i, j);
end;

% plot dipole for one specific dataset
figure;
for DATASET = [6 8 14 2 9]
    switch DATASET
     case 1 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/km81/Sternberg/');
     case 2 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/jo74/Sternberg/');
     case 3 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ds76/Sternberg/');
     case 4 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/cj82/Sternberg/');
     case 5 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ap82/Sternberg/');
     case 6 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ke70/Sternberg/');
     case 7 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/tp62/Sternberg/');
     case 8 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/cz84/Sternberg/');
     case 9 , EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/gm84/Sternberg/');
     case 10, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/gv84/Sternberg/');
     case 11, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/nf68/Sternberg/');
     case 12, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ds80/Sternberg/');
     case 13, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/kb77/Sternberg/');
     case 14, EEG = pop_loadset( 'tmp4ica.set', '/data/common1/stern/eeg/ts79/Sternberg/');
    end;
    filename2 = sprintf('randomdip_%d_data.mat',DATASET);
    load('-mat', filename2);

    indices = round(rand(1, 90)*(EEG.trials*EEG.pnts-1))+1;
    
    count = 1;
    for index = indices
        subplot(10,10,count);
        topoplot(double(EEG.data(:,index)), EEG.chanlocs, 'electrodes', 'off');
        count=count+1;
    end;
    textsc([ 'Dataset ' int2str(DATASET)], 'title');
    set(gcf, 'paperpositionmode', 'auto', 'position', [520   243   874   867]);
    eval(sprintf('print -depsc dataset_dipoles%d.eps', DATASET));
    clf
end;
