%% PCA via SVD of covariance
%% Simulate data with covariance structure
% Simulation parameters
N = 10000; % time points
M =    20; % channels

% Time vector (radian units)
t = linspace(0,6*pi,N);

% Relationship across channels (imposing covariance)
chanrel = sin(linspace(0,2*pi,M))';

% Generate data
data = bsxfun(@times,repmat( sin(t),M,1 ),chanrel) + randn(M,N);

%% PCA via eigendecomposition of covariance matrix
% First mean-center (mean over time!)
dataM = data - mean(data,2);

% All pairwise dot products as the matrix times its transpose
covmat = dataM*dataM'/(N-1);

% Sort
[Veig,Deig] = eig( covmat );
[Deig,sidx] = sort( diag(Deig),'descend' );
Veig = Veig(:,sidx);

% Convert eigenvalues to %
Deig = 100*Deig/sum(Deig);

% Time series of top component
eig_ts = Veig(:,1)'*data;

%% PCA via SVD of data matrix
% SVD
[U,S,V] = svd( dataM,'econ' );

% Convert singular values to %
S = S.^2; % makes it comparable to eigenvalues
s = 100*diag(S)/sum(diag(S));

%% Show the results
figure(1), clf

% Plot eigenvalue/singular-value spectrum
subplot(221), hold on
plot(Deig,'bs-','markerfacecolor','w','markersize',10,'linew',2)
plot(s,'ro-','markerfacecolor','w','markersize',5)
xlabel('Component number'), ylabel('\lambda or \sigma')
legend({'eig';'svd'})
title('Eigenspectrum')

% Show eigenvector/singular value
subplot(222), hold on
plot(Veig(:,1),'bs-','markerfacecolor','w','markersize',10,'linew',2)
plot(-U(:,1),'ro-','markerfacecolor','w','markersize',5)
xlabel('Vector component'), ylabel('Weight')
title('Component weights')

% Time series
subplot(212)
timevec = (0:N-1)/N;
plot(timevec,eig_ts, timevec,-V(:,1)*sqrt(S(1)),'.')
xlabel('Time (norm.)')
title('Component time series')
legend({'eig';'svd'})
zoom on

%% end.