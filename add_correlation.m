dat = 2;
EEG = geticadata(dat, 'Ext. Infomax');
%EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data(:,:);

lng = size(EEG.icaact,1);
[MI,vMI] = minfo(EEG.icaact(:,:));
PMI = (sum(MI(:)) - sum(diag(MI))) / (lng*(lng-1));

return;

pmi2 = load('-mat', 'pmi_save_data.mat');
tmppmi = PMI/pmi2.PMI(dat)*100;

% add correlation to components 14 and 35 (two mu rythms)
tmpdata = EEG.icaact([14 35], :);
pmival  = minfo(tmpdata);

% -----------
% DATA PMI + 40 Hz
dat = 2;
EEG = geticadata(dat, 'Amica');
lng = EEG.nbchan;
[MI2,vMI2] = minfo(EEG.data(:,:));
PMI = (sum(MI(:)) - sum(diag(MI))) / (lng*(lng-1));

EEG.data(:,150:168,:) = EEG.data(:,150:168,:) + repmat(sin(linspace(0, 6*pi, 19))*std(EEG.data(:)), [EEG.nbchan 1 EEG.trials]);
lng = EEG.nbchan;
[MI2,vMI2] = minfo(EEG.data(:,:));
PMI2 = (sum(MI2(:)) - sum(diag(MI2))) / (lng*(lng-1));

% -----------
% COMPONENT PMI + 40 Hz
dat = 2;
EEG = geticadata(dat, 'Amica');
lng = size(EEG.icaact,1);
[MI,vMI] = minfo(EEG.icaact(:,:));
PMI = (sum(MI(:)) - sum(diag(MI))) / (lng*(lng-1));

tmpica = EEG.icaact(1:2,:,:);
EEG.icaact(1:2,150:168,:) = EEG.icaact(1:2,150:168,:) + repmat(sin(linspace(0, 6*pi, 19))*std(tmpica(:)), [2 1 EEG.trials]);
[MI2,vMI2] = minfo(EEG.icaact(:,:));
PMI2 = (sum(MI2(:)) - sum(diag(MI2))) / (lng*(lng-1));

