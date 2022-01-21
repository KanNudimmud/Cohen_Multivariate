%% Overfitting and inferential statistics 
%  What is overfitting and why is it inappropriate?
%% Polynomial example
N = 10;

% Generate data as linear trend plus noise
data1 = linspace(-2,2,N) + randn(1,N);
data2 = linspace(-2,2,N) + randn(1,N);

% Polynomial fit, based on data1, with 1 or N parameters
p1 = polyval( polyfit(1:N,data1,1),1:N );
p2 = polyval( polyfit(1:N,data1,N),1:N );

% Plot dataset 1 and its fit
figure(1), clf

subplot(121), hold on
plot(data1,'rs','markerfacecolor','k','markersize',15)
plot(p1,'bo-','linew',2,'markersize',15,'markerfacecolor','w')
plot(p2,'kp-','linew',2,'markersize',15,'markerfacecolor','w')
title('Data1')
legend({'Data';'Low-D fit';'High-D fit'})
axis square

subplot(122), hold on
plot(data2,'rs','markerfacecolor','k','markersize',15)
plot(p1,'bo-','linew',2,'markersize',15,'markerfacecolor','w')
plot(p2+.1,'kp-','linew',2,'markersize',15,'markerfacecolor','w')
title('Data2')
axis square

%% For "EEG" data
% Mat file containing EEG, leadfield and channel locations
load emptyEEG
EEG.srate = 500;

EEG.pnts  = 1000; % time points
EEG.times = (0:EEG.pnts-1)/EEG.srate;
EEG.data  = squeeze(lf.Gain(:,1,:))*randn(size(lf.Gain,3),EEG.pnts);

%% Covariance matrices and GED
% Time period 1
tmpdat  = EEG.data(:,1:round(EEG.pnts/2));
tmpdat  = bsxfun(@minus,tmpdat,mean(tmpdat,2));
covmat1 = tmpdat*tmpdat' / EEG.pnts/2;

% Time period 2
tmpdat  = EEG.data(:,round(EEG.pnts/2)+1:end);
tmpdat  = bsxfun(@minus,tmpdat,mean(tmpdat,2));
covmat2 = tmpdat*tmpdat' / EEG.pnts/2;

% GED
[W,L] = eig( covmat1,covmat2 );
[l,sidx] = sort( diag(L),'descend' );
W = W(:,sidx);

%% Plotting
% Eigenvalues
figure(2), clf
subplot(211)
plot( l,'ks-','linew',2,'markerfacecolor','w','markersize',15 )
set(gca,'xlim',[0 EEG.nbchan+1])
xlabel('Component number'), ylabel('\lambda')
title('Eigenspectrum')

% Topographical maps
clim = [-1 1]*100;
for i=1:3
    subplot(2,3,3+i)
    topoplotIndie(W(:,i)'*covmat1,EEG.chanlocs);
    set(gca,'clim',clim)
    title([ 'Component ' num2str(i) ])
end
colormap jet

%% end.