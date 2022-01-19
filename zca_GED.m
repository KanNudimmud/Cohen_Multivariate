%% Source separation with GED
%  ZCA+two-stage separation in real data
%% Load data
load sampleEEGdata.mat
EEG.data = double( EEG.data );

%% ZCA
tidx = dsearchn(EEG.times',[-1000 0]');

% Covariance of pre-stim activity
tmpdat = reshape( EEG.data(:,tidx(1):tidx(2),:), [ EEG.nbchan (diff(tidx)+1)*EEG.trials ]);
tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
covN   = tmpdat*tmpdat'/size(tmpdat,2);

% Eig and convert to %change
[V,D] = eig(covN);
d = diag(D);
d = 100*d/sum(d);

% Threshold and replace data with whitened data
r = d>.1;
EEG.data = reshape( V(:,r)*D(r,r)^(-1/2)*V(:,r)'*reshape(EEG.data,EEG.nbchan,[]) ,size(EEG.data) );

%% Stage 1: Compression
% Covariance from all broadband data
covA = zeros(EEG.nbchan);
for triali=1:EEG.trials
    tmp = detrend( EEG.data(:,:,triali)' )';
    covA = covA + tmp*tmp'/EEG.pnts;
end

% PCA and sort
[V,D]  = eig(covA/EEG.trials);
[D,sx] = sort(diag(D),'descend');

% Sort eigvecs
V = V(:,sx);

% Convert evals to %change
D = 100*D/sum(D);

% Which components to keep (threshold of % total variance accounted for)
comps2keep = D>1;
nComps = sum(comps2keep);

% Compute subspaced data
EEG.subdata = reshape( (reshape(EEG.data,EEG.nbchan,[])' * V(:,comps2keep))' ,[nComps EEG.pnts EEG.trials]);

%% Stage 2
% Filter PC components (compressed data)
EEG.fdata = filterFGx(EEG.subdata,EEG.srate,11,4);

% Covariance matrices
[covT,covB] = deal( zeros(nComps) );
for i=1:EEG.trials
    tmp = detrend(squeeze(EEG.subdata(:,:,i))')';
    covB = covB + tmp*tmp'/EEG.pnts;
    
    tmp = detrend(squeeze(EEG.fdata(:,:,i))')';
    covT = covT + tmp*tmp'/EEG.pnts;
end

% GED the shit out that mug, yo.
[evecs,evals] = eig(covT,covB);
[evals,sidx]  = sort( diag(evals),'descend' );

% Time series from only the PCA components
ts = reshape( reshape(EEG.subdata,nComps,[])'*evecs(:,sidx(1)) ,EEG.pnts,EEG.trials);

% Note the additional transformation using PCA eigenvectors
topo = evecs(:,sidx(1))'*covT*V(:,1:nComps)';
[~,mv] = max(abs(topo));
if topo(mv)<0
    topo = -topo;
    ts = -ts;
end

%% Plot the results
figure(1), clf

% Topomap
subplot(221)
topoplotIndie(topo,EEG.chanlocs,'numcontour',0);
title('Filter forward model of top component')

% Eigspect
subplot(222)
plot(evals,'s-','markerfacecolor','w','linew',2,'markersize',10)
xlabel('Component number'), ylabel('\lambda')
title('Eigenspectrum')
box off

% Time series
subplot(212)
plot(EEG.times,mean(ts,2),'linew',3)
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)')
title('Component time series')
box off

%% end.