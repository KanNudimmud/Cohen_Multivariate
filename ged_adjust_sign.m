%% Source separation with GED
% Adjusting component signs
%% Preliminaries
% Mat file containing EEG, leadfield and channel locations
load emptyEEG
EEG.srate = 500;

EEG.trials = 1;
EEG.pnts   = 1000; % this is hard-coded below, so don't change this :/
EEG.times  = (0:EEG.pnts-1)/EEG.srate;
EEG.data   = zeros(EEG.nbchan,EEG.pnts,EEG.trials);

%% Ground-truth dipole projection
diploc = 109;

figure(1), clf
clim = [-45 45];
subplot(231)
topoplotIndie(-lf.Gain(:,1,diploc), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Ground truth map')
xlabel('Time (s)'), ylabel('Activity (a.u.)')

colormap jet

%% Insert activity waveforms into dipole data
% Now we generate random data in brain dipoles.
% create 1000 time points of random data in brain dipoles
% (note: the '1' before randn controls the amount of noise)
dipole_data = 1*randn(length(lf.Gain),1000);

% Add signal to second half of dataset
dipole_data(diploc,501:end) = -15*sin(2*pi*10*(0:499)/EEG.srate);

% Normalize dipoles (simplifies the model)
lf.GainN = bsxfun(@times,squeeze(lf.Gain(:,1,:)),lf.GridOrient(:,1)') + bsxfun(@times,squeeze(lf.Gain(:,2,:)),lf.GridOrient(:,2)') + bsxfun(@times,squeeze(lf.Gain(:,3,:)),lf.GridOrient(:,3)');

% Project dipole data to scalp electrodes
EEG.data = lf.GainN*dipole_data;

% Plot the time series
subplot(212), hold on
plot(EEG.times,dipole_data(diploc,:)/2,'linew',2)

%% Create covariance matrices
% Compute covariance matrix R is first half of data
tmpd = EEG.data(:,1:500);
tmpd = bsxfun(@minus,tmpd,mean(tmpd,2));
covR = tmpd*tmpd'/500;

% Compute covariance matrix S is second half of data
tmpd = EEG.data(:,501:end);
tmpd = bsxfun(@minus,tmpd,mean(tmpd,2));
covS = tmpd*tmpd'/500;

% GED and sort solutions
[evecs,evals] = eig(covS,covR);
[evals,sidx]  = sort(diag(evals),'descend');
evecs = evecs(:,sidx);

%% Option 1: ignore the sign
% Compute filter forward model and component time series
map   = evecs(:,1)'*covS;
cdat1 = evecs(:,1)'*EEG.data;

% Visualize
subplot(232)
topoplotIndie(map, EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Option 1 (ignore)')

subplot(212)
plot(EEG.times,cdat1,'linew',2)

%% Option 2: force positive max. value
% Find max magnitude and its sign
[~,idx] = max(abs(map));
signOfMax = sign(map(idx));

% Sign-flip the eigenvector
w = evecs(:,1)*signOfMax;

% Recompute the map and component time series
map = w'*covS;
cdat2 = w'*EEG.data;

% Visualize
subplot(233)
topoplotIndie(map, EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','off','shading','interp');
title('Option 2 (max-pos)')

subplot(212)
h = plot(EEG.times,cdat2,'linew',2); % with handle for 2b below

legend({'Dipole data','Option 1','Option 2'})
set(gca,'xlim',[1 2])

%% Option 2b: force sign of a particular electrode
chan2bPos = 'Fz';

% Find max magnitude and its sign
map = evecs(:,1)'*covS;
targetSign = sign(map( strcmpi({EEG.chanlocs.labels},chan2bPos) ));

% Sign-flip the eigenvector
w = evecs(:,1)*targetSign;

% Recompute the map and component time series
map = w'*covS;
cdat2 = w'*EEG.data;

%% Visualize
subplot(233)
topoplotIndie(map, EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','labels','shading','interp');
title('Option 2b (chan-pos)')

subplot(212)
set(h,'ydata',cdat2)

legend({'Dipole data','Option 1','Option 2b'})
set(gca,'xlim',[1 2])

%% end