%% PCA on multi-trial data
%% Load data
load sampleEEGdata.mat
EEG.data = double( EEG.data );

%% Erp ("phase-locked") covariance matrix
% Time window to use
tidx = dsearchn(EEG.times',[0 800]');

% Phase-locked covariance
erp  = mean(EEG.data(:,tidx(1):tidx(2),:),3);
data = erp - mean(erp,2);
cov_phaselocked = data*data'/(EEG.pnts-1);

figure(1), clf
subplot(121)
imagesc(cov_phaselocked)
axis square
set(gca,'clim',[-1 1]*2)
xlabel('Channels'), ylabel('Channels')
colorbar
title('Covariance of average')

%% Single-trial ("total") covariance matrix
cov_total = zeros(EEG.nbchan);

for triali=1:EEG.trials
    % Extract data and compute mini-cov
    minicov = cov( EEG.data(:,tidx(1):tidx(2),triali)' );
    
    % Sum into covariance
    cov_total = cov_total + minicov;
end

% Scale by N
cov_total = cov_total/triali;

% Visualize
subplot(122)
imagesc(cov_total)
axis square
set(gca,'clim',[-1 1]*100)
xlabel('Channels'), ylabel('Channels')
colorbar
title('Average of covariances')

%% Compute and visualize phase-locked PCA
% PCA of trial-average
[V_PL,L_PL] = eig( cov_phaselocked );
[L_PL,sidx] = sort(diag(L_PL),'descend');
V_PL = V_PL(:,sidx);
L_PL = 100*L_PL/sum(L_PL);

% Principal components time series
% NOTE: applying filter to the entire time series!
compts_PL = V_PL(:,1)' * mean(EEG.data,3);

% Visualize
figure(2), clf
subplot(321)
plot(L_PL,'ks-','markerfacecolor','w','markersize',8)
set(gca,'xlim',[.5 20])
xlabel('Component number')
ylabel('Pct variance')
title('Scree plot')

subplot(322)
plot(cumsum(L_PL),'ks-','markerfacecolor','w','markersize',8)
set(gca,'xlim',[.5 20])
xlabel('Component number')
ylabel('Cumulative sum variance')

% Topoplot
subplot(323)
topoplotIndie(V_PL(:,1),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','off');
title('Phase-locked c_1')
colormap jet

% Time course
subplot(313)
plot(EEG.times,compts_PL,'linew',2)
set(gca,'xlim',[-200 1300])
xlabel('Time (ms)')
ylabel('PC activity')

%% Repeat for total PCA
% PCA of trial-average
[V_TT,L_TT] = eig( cov_total );
[L_TT,sidx] = sort(diag(L_TT),'descend');
V_TT = V_TT(:,sidx);
L_TT = 100*L_TT/sum(L_TT);

% Principal components time series
% data are 3D, so we reshape to 2D, project, then reshape back
data2d    = reshape(EEG.data,EEG.nbchan,[]);
comptmp   = V_TT(:,1)'*data2d;
compts_TT = reshape(comptmp,[ EEG.pnts EEG.trials ]);

% Then compute the trial average
compTT_erp = mean(compts_TT,2);

% Visualize
subplot(321), hold on
plot(L_TT,'ro-','markerfacecolor','w','markersize',8)
legend({'Phase-locked','Total'})

subplot(322), hold on
plot(cumsum(L_TT),'ro-','markerfacecolor','w','markersize',8)
legend({'Phase-locked','Total'})

% Topoplot
subplot(324)
topoplotIndie(V_TT(:,1),EEG.chanlocs,'numcontour',0,'shading','interp','electrodes','off');
title('Total c_1')

% Time course
subplot(313), hold on
plot(EEG.times,compTT_erp,'linew',2)
legend({'Phase-locked','Total'})

%% end