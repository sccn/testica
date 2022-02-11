algorithms = {  'Amica' 'Ext. Infomax' 'Pearson' 'Infomax' ...
                 'SHIBBS'   'FastICA'    'JADE'    'TICA'    'JADE opt.'    'JADE-TD' ...
                 'FOBI'    'SOBIRO'    'EVD24'    'EVD'    'SOBI'    'icaMS' ...
                'AMUSE'  'PCA' 'sphere' };
nchans = 71;
h = zeros(nchans,length(algorithms),14);

clear mir
for algo=1:length(algorithms)
   for dat=1:14
      EEG = geticadata(dat, algorithms{algo});
      W = EEG.icaweights;
      S = EEG.icasphere;
      WS = W*S;
      
      %load('-mat',['ica' int2str(dat) '_72_14_simbec.mat']);
      %WS = W;
      %DATASET = dat;
      %processdat;

      if algo == 1
          %h0(:,dat) = getent(reshape(EEG.data,nchans,EEG.pnts*EEG.trials));
          h0(:,dat) = getent2(reshape(EEG.data,nchans,EEG.pnts*EEG.trials));
      end
      s = WS * reshape(EEG.data,nchans,EEG.pnts*EEG.trials); 
      h(:,algo,dat) = getent2(s);
      %h(:,algo,dat) = getent(s);
      %mir(algo,dat) = sum(h0(:,dat)) - sum(h(:,algo,dat)) + log(abs(det(WS)))
      mir(algo,dat) = sum(h0(:,dat)) - sum(h(:,algo,dat)) + sum(log(abs(eig(WS))))
   end;
end;

% save results
% WARNING: THIS IS GOING TO ERASE THE EXISTING FILE
copyfile('mir_new.mat', 'mir_new_old.mat');
save('-mat', 'mir_new.mat', 'mir', 'algorithms');
