%% Source separation with GED
% Two components in simulated EEG data
%% Preliminaries
% Mat file containing EEG, leadfield and channel locations
load emptyEEG
EEG.srate = 500;

EEG.trials = 200;  % total, 1/2 per condition
EEG.pnts   = 1000; % time points per trial
EEG.times  = (0:EEG.pnts-1)/EEG.srate;
EEG.data   = zeros(EEG.nbchan,EEG.pnts,EEG.trials);

%% Dipole Locations
dipoleLoc1 = 109;
dipoleLoc2 = 135;

figure(1), clf
clim = [-45 45];
subplot(331)
topoplotIndie(-lf.Gain(:,1,dipoleLoc1), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','numbers','shading','interp');
title('Simulation dipole 1')

subplot(332)
topoplotIndie(-lf.Gain(:,1,dipoleLoc2), EEG.chanlocs,'maplimits',clim,'numcontour',0,'electrodes','numbers','shading','interp');
title('Simulation dipole 2')

colormap jet

%% Insert Activity Waveforms into Dipole Data
% Frequencies of the two dipoles
freq1 = 15;
freq2 = 10;

% Time point of "stimulus" onset
tidx = dsearchn(EEG.times',mean(EEG.times));

% The "innards" of the sine function
omega1 = 2*pi*freq1*EEG.times(tidx:end);
omega2 = 2*pi*freq2*EEG.times(tidx:end);

% Loop over trials
for ti=1:EEG.trials
    % Source waveforms (sine waves with random phase)
    swave1 = sin( omega1 + rand*2*pi );
    swave2 = sin( omega2 + rand*2*pi );
    
    dipole_data = randn(size(lf.Gain,3),EEG.pnts)/5;
    dipole_data(dipoleLoc1,tidx:end) = dipole_data(dipoleLoc1,tidx:end) + swave1;
    dipole_data(dipoleLoc2,tidx:end) = dipole_data(dipoleLoc2,tidx:end) + swave2;
    
    % Project to scalp
    EEG.data(:,:,ti) = squeeze(lf.Gain(:,1,:))*dipole_data;
end

%% GED for Spatial Filter
[covPre,covPst] = deal( zeros(EEG.nbchan) );

% Covariance matrices per trial
for ti=1:EEG.trials
    % "prestim" covariance
    tdat   = squeeze(EEG.data(:,1:tidx,ti));
    tdat   = tdat - mean(tdat,2);
    covPre = covPre + (tdat*tdat') / 500;
    
    % "post-stim" covariance
    tdat   = squeeze(EEG.data(:,tidx:end,ti));
    tdat   = tdat - mean(tdat,2);
    covPst = covPst + (tdat*tdat') / 500;
end

covPre = covPre./ti;
covPst = covPst./ti;

% GED
[evecs,evals] = eig(covPst,covPre);
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);

%%% Compute filter forward models and flip sign
% Component 1:
maps(:,1) = evecs(:,1)'*covPst; % get component
[~,idx]   = max(abs(maps(:,1)));  % find max magnitude
maps(:,1) = maps(:,1)*sign(maps(idx,1)); % possible sign flip

% Component 2:
maps(:,2) = evecs(:,2)'*covPst; % get component
[~,idx]   = max(abs(maps(:,2)));  % find max magnitude
maps(:,2) = maps(:,2)*sign(maps(idx,2)); % possible sign flip

%%% Compute component time series (projections)
cdat = evecs(:,1:2)'*reshape(EEG.data,EEG.nbchan,[]);
cdat = reshape( cdat, [ 2 EEG.pnts EEG.trials ]);

%% Show Topographical Maps and Eigenspectrum
for i=1:2
    subplot(3,3,3+i)
    topoplotIndie(maps(:,i),EEG.chanlocs,'numcontour',0,'electrodes','off');
    title([ 'Component ' num2str(i) ])
end

subplot(336), plot(evals,'s-','linew',2,'markersize',10,'markerfacecolor','k')
set(gca,'xlim',[0 15])
ylabel('\lambda'), title('Eigenvalues of decomposition'), axis square

%% Time-Frequency Analysis on Components
% Frequencies in Hz
frex = linspace(2,20,20);

% Convenient to have component time series data as 2D
comp2d = reshape(cdat,2,[]);
% comp2d = reshape(EEG.data([31 57],:,:),2,[]);

% Initialize time-frequency matrix
ctf = zeros(2,length(frex),EEG.pnts);

% Loop over frequencies
for fi=1:length(frex)
    % Filter data for both components at this frequency
    filtdat = filterFGx(comp2d,EEG.srate,frex(fi),4);
    
    % Loop over components
    for compi=1:2
        % Compute power time series as envelope of Hilbert transform
        as = reshape(hilbert(filtdat(compi,:)) ,EEG.pnts,EEG.trials);
        
        % TF power is trial-average power
        ctf(compi,fi,:) = mean( abs(as).^2 ,2);
    end
end

%%  Plotting
for compi=1:2
    subplot(3,3,6+compi)
    
    contourf(EEG.times,frex,squeeze(ctf(compi,:,:)),40,'linecolor','none')
    
    axis square
    title([ 'Component ' num2str(compi) ])
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
end

%% end