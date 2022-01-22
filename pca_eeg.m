%% PCA of simulated EEG data
%% Simulate data
% mat file containing EEG, leadfield and channel locations
load emptyEEG

% Index of dipole to simulate activity in
diploc = 109;

% Plot brain dipoles
figure(1), clf, subplot(121)
plot3(lf.GridLoc(:,1), lf.GridLoc(:,2), lf.GridLoc(:,3), 'bo','markerfacecolor','y')
hold on
plot3(lf.GridLoc(diploc,1), lf.GridLoc(diploc,2), lf.GridLoc(diploc,3), 'rs','markerfacecolor','k','markersize',10)
rotate3d on, axis square
title('Brain dipole locations')

% Each dipole can be projected onto the scalp using the forward model. 
% The code below shows this projection from one dipole.
subplot(122)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'electrodes','numbers','shading','interp');
title('Signal dipole projection')

% Number of time points and time vector
N = 1000;
EEG.times = (0:N-1)/EEG.srate;

% Random data in all brain dipoles
dipole_data = randn(length(lf.Gain),N);

% Add signal to one dipole
dipole_data(diploc,:) = 15*sin(2*pi*10*EEG.times);

% Project data from all dipoles to scalp electrodes
EEG.data = squeeze(lf.Gain(:,1,:))*dipole_data;

%% Compute PCA
% Mean-center data
data = EEG.data - mean(EEG.data,2);

% Covariance matrix
covd = data*data'/size(EEG.data,2);

% Eigendecomposition
[evecs,evals] = eig( covd );

% Sort according to eigenvalues
[evals,sidx] = sort(diag(evals),'descend');
evecs = evecs(:,sidx);
evals = 100*evals/sum(evals);

% Principal component time series
pc_timeseries = evecs(:,1)'*EEG.data;

%% Plot results and compare with ground truth
figure(2), clf

% Eigenspectrum
subplot(231)
plot(evals(1:20),'ko-','markerfacecolor','w','linew',2)
title('Eigenspectrum'), axis square
ylabel('Percent variance'), xlabel('Component number')

% Topographical map of first eigenvector
subplot(232)
topoplotIndie(evecs(:,1),EEG.chanlocs,'numcontour',0,'shading','interp');
title('PC topomap')

% Topographical map of dipole (ground truth)
subplot(233)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'numcontour',0,'shading','interp');
title('Ground truth topomap')

% Plot time series
subplot(212)
plot(EEG.times,pc_timeseries, ...
     EEG.times,EEG.data(31,:), ...
     EEG.times,-dipole_data(diploc,:)*100 ,'linew',2 )
legend({'PC';'Chan';'Dipole'})
xlabel('Time (s)')

%% end.