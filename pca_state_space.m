%% State-space representation via PCA
%% Load data
load sampleEEGdata.mat
EEG.data = double( EEG.data );

%% Compute PCA
% Average data, mean-center, covariance
erp  = mean(EEG.data,3);
data = erp - mean(erp,2);
covd = data*data'/(EEG.pnts-1);

% Eig and sort
[evecs,evals] = eig( covd );
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);
evals = 100*evals/sum(evals);

% Principal components time series
pc_timeseries = evecs(:,1:2)'*erp;

%% Plot PC results in time-voltage space
figure(1), clf

% Eigenspectrum
subplot(231)
plot(evals(1:20),'ko-','markerfacecolor','w','linew',2)
title('Eigenspectrum'), axis square
ylabel('Percent variance'), xlabel('Component number')

% Topographical map of first eigenvector
subplot(232)
topoplotIndie(evecs(:,1),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','numbers');
title('PC1 topomap')

subplot(233)
topoplotIndie(evecs(:,2),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','numbers');
title('PC2 topomap')

% Plot time series
subplot(212)
plot(EEG.times,pc_timeseries,'linew',2 )
legend({'PC1';'PC2'})
xlabel('Time (s)')

%% Plot state-space representation
% Time windows for plotting in ms
timewin = [0 250; 200 600; 600 800];

% Smooth the PC time series
pcts = pc_timeseries;
k = 5; % smoothing kernel
for ti=k+1:EEG.pnts-k-1
    pcts(1,ti) = mean(pc_timeseries(1,ti-k:ti+k));
    pcts(2,ti) = mean(pc_timeseries(2,ti-k:ti+k));
end

% Plot the data
figure(2), clf, hold on
leg = cell(1,size(timewin,1));
for twini=1:size(timewin,1)
    % Convert to indices
    tidx = dsearchn(EEG.times',timewin(twini,:)');
    
    % Plot
    plot(pcts(1,tidx(1):tidx(2)),pcts(2,tidx(1):tidx(2)),'linew',3)
    
    % Define legend
    leg{twini} = [ num2str(timewin(twini,1)) '-' num2str(timewin(twini,2)) ' ms' ];
end

xlabel('PC 1'), ylabel('PC 2')
legend(leg)

%% end.