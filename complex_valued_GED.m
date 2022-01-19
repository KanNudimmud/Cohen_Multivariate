%% Source separation with GED
% Complex-valued solutions
%% Load Data
load lowrankEEG.mat
EEG.data = double(EEG.data);

%% Rank of Matrix
r = rank(EEG.data(:,:,1));
disp("rank = " + r)

%% Get Covariance Matrices
% Filter data in alpha
EEG.fdata = filterFGx(EEG.data,EEG.srate,11,4);

% Initialize to zeros
[covS,covR] = deal( zeros(EEG.nbchan) );

for i=1:EEG.trials
    % Detrend instead of demean
    tmp = detrend(squeeze(EEG.data(:,:,i))')';
    covR = covR + tmp*tmp'/EEG.pnts;
    
    tmp = detrend(squeeze(EEG.fdata(:,:,i))')';
    covS = covS + tmp*tmp'/EEG.pnts;
end 

%% Regularized GED
% Regularization parameter
regu_gam = .0;

% Apply regularization
R = (1-regu_gam)*covR + regu_gam*mean(eig(covR))*eye(EEG.nbchan);

% Report ranks of matrices
clc
disp("rank(S)  = " + rank(covS))
disp("rank(R)  = " + rank(covR))
disp("rank(Rr) = " + rank(R))

% GED etc
[evecs,evals] = eig(covS,R);
[evals,sidx]  = sort( diag(evals),'descend' );
evecs = evecs(:,sidx);
ts = reshape( reshape(EEG.data,EEG.nbchan,[])'*evecs(:,1) ,EEG.pnts,EEG.trials);

% Adjust sign
topo = evecs(:,1)'*covS;
[~,mv] = max(abs(topo));
if topo(mv)<0
    topo = -topo;
    ts = -ts;
end

%% Plotting
figure(1), clf
subplot(211)
topoplotIndie(topo,EEG.chanlocs);

subplot(212)
plot(EEG.times,mean(ts,2),'k','linew',3)
set(gca,'xlim',[-200 1000])
xlabel('Time (ms)')

%% end