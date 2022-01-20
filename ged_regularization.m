%% Source separation with GED
% Effects of regularization in real data
%% 
% The transformation matrix T
T = eye(3);
T(end) = 0;

% Create the spherical coordinates
[x,y,z] = sphere;
M = [x(:) y(:) z(:)]';

%%% Setup the figure
figure(1), clf
h = surf(x,y,z);
shading interp
axis([-1 1 -1 1 -1 1])
axis square, rotate3d on
view([-70 15])

% Loop over gamma values
for gami=linspace(0,.5,20)
    Tt = T*(1-gami) + gami*mean(eig(T))*eye(3);
    
    % Impose the transformation matrix onto the sphere coordinates
    TM = Tt*M;
    
    % Reshape back to new coordinatess
    xp = reshape(TM(1,:),[21 21]); % note: 21 is default resolution for sphere function
    yp = reshape(TM(2,:),[21 21]);
    zp = reshape(TM(3,:),[21 21]);
    
    % Update the plot
    set(h,'XData',xp,'YData',yp,'ZData',zp)
    title([ 'Regularization \gamma = ' num2str(gami) ])
    pause(.3)
end

%% Now for real data...

load sampleEEGdata.mat

% Regularization parameter
regu_gam = .0;

% Filter data in alpha
EEG.fdata = filterFGx(EEG.data,EEG.srate,11,4);

% Get covariances from each trial
[covF,covB] = deal( zeros(EEG.nbchan) );

for i=1:EEG.trials
    tmp  = detrend(squeeze(EEG.data(:,:,i))')';
    covB = covB + tmp*tmp'/EEG.pnts;
    
    tmp  = detrend(squeeze(EEG.fdata(:,:,i))')';
    covF = covF + tmp*tmp'/EEG.pnts;
end 

covB = covB/i;
covF = covF/i;

% Apply regularization
Rregu = (1-regu_gam)*covB + regu_gam*mean(eig(covB))*eye(EEG.nbchan);

% GED
[evecs,evals] = eig(covF,Rregu);
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);
ts = reshape( reshape(EEG.data,EEG.nbchan,[])'*evecs(:,1) ,EEG.pnts,EEG.trials);

% Topomaps and time series
topo = evecs(:,1)'*covF;
[~,mv] = max(abs(topo));
if topo(mv)<0
    topo = -topo;
    ts = -ts;
end

figure(2)
subplot(311)
topoplotIndie(topo,EEG.chanlocs);

subplot(312)
plot(EEG.times,mean(ts,2),'k','linew',3)
set(gca,'xlim',[-200 1000])

%%% Show both covariance matrices
subplot(3,3,7)
imagesc(covB)
axis square
set(gca,'xtick',[],'ytick',[],'clim',[-1 1]*150)
title('R matrix')
colorbar

subplot(3,3,8)
imagesc(Rregu)
set(gca,'xtick',[],'ytick',[],'clim',[-1 1]*150)
axis square
title([ 'Regularized (\gamma=' num2str(100*regu_gam) '%)' ])
colorbar

subplot(3,3,9)
imagesc(covB-Rregu)
set(gca,'xtick',[],'ytick',[],'clim',[-1 1]*10)
axis square
title('R - Rregu')
colorbar
colormap jet

%% end.