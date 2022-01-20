%% Two-stage source separation in real EEG data
% Load the data
load sampleEEGdata.mat
EEG.data = double(EEG.data);

%% Compute PCA on ERP
% Initialize the covariance matrix
covA = zeros(EEG.nbchan);

% Compute the covariance matrix for each trial
for triali=1:EEG.trials
    tmp  = EEG.data(:,:,triali);
    covA = covA + tmp * tmp' / (EEG.pnts-1);
end

% Take mean of covariance matrix
covA = covA / EEG.trials;

% Calculate eigenvectors and values
[evecs,evals] = eig(covA);

% Sort eigenvalues
[evals,sidx] = sort(diag(evals),'descend');
evecs        = evecs(:,sidx);

% Convert eigenvalues to percentage
evals = 100 * evals / sum(evals);

% Components accounting for more than 1%
comps2keep = evals > 1;
nComps     = sum(comps2keep);

% Principal component time series with the component accounting for
pcts = zeros(nComps,EEG.pnts,EEG.trials);
for triali=1:EEG.trials
    pcts(:,:,triali) = evecs(:,comps2keep)' * EEG.data(:,:,triali);
end

%% GED
% Apply a narrow band-pass filter centered at 11 Hz. with a FWHM of 4 Hz.
fdata = filterFGx(pcts,EEG.srate,11,4);

% For each trial, compute the S matrix
fcov = zeros(nComps);
for triali=1:EEG.trials
    tmp  = fdata(:,:,triali);
    fcov = fcov + tmp * tmp' / (EEG.pnts-1);
end

% Get the mean covariance matrix
fcov = fcov/EEG.trials;

% For each trial, compute the R matrix
bcov = zeros(nComps);
for triali=1:EEG.trials
    tmp  = pcts(:,:,triali);
    bcov = bcov + tmp * tmp' / (EEG.pnts-1);
end

% get the mean covariance matrix
bcov = bcov / EEG.trials;

% Get eigenvectors and eigenvalues with GED
[g_evecs,g_evals] = eig(fcov,bcov);
[g_evals,sidx]    = sort(diag(g_evals),'descend');
g_evecs           = g_evecs(:,sidx);

% Component time series from only the PCA components
gts = zeros(EEG.pnts,EEG.trials);
for triali= 1:EEG.trials
    gts(:,triali) = g_evecs(:,1)' * pcts(:,:,triali);
end

% Project map through the biggest eigenvector
topomap = g_evecs(:,1)' * fcov * evecs(:,comps2keep)';

% Adjust sign
[~,maxi] = max(abs(topomap));
gts      = gts * sign(topomap(maxi));
topomap  = topomap * sign(topomap(maxi));

% Plot
figure(2), clf
subplot(211)
topoplotIndie(topomap,EEG.chanlocs,'numcontour',0,'electrodes','off');

subplot(212)
plot(EEG.times,mean(gts,2),'k','linew',3);
xlabel('Time (ms)'), ylabel('Activity')
set(gca,'xlim',[-200,1000],'fontsize',18)

%% end.