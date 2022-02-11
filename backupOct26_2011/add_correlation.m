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

