%% How to perform a principal components analysis
%% Simulate data with covariance structure
% Simulation parameters
N = 1000;     % time points
M =   20;     % channels
nTrials = 50; % number of trials

% Time vector (radian units)
t = linspace(0,6*pi,N);

% Relationship across channels (imposing covariance)
chanrel = sin(linspace(0,2*pi,M))';

data = bsxfun(@times,repmat( sin(t),M,1 ),chanrel) + randn(M,N);

% Step 1: mean-center and compute covariance matrix
data = bsxfun(@minus,data,mean(data,2));
covmat = data*data'/(N-1);

% Step 2: eigendecomposition
[evecs,evals] = eig( covmat );

% Step 3: sort vectors by values
[evals,soidx] = sort( diag(evals),'descend' );
evecs = evecs(:,soidx);

% Step 4: compute component time series
r = 2; % two components
comp_time_series = evecs(:,1:r)'*data;

% Step 5: convert eigenvalues to percent change
evals = 100*evals./sum(evals);

% Step 6: visualize and interpret the results
figure(1), clf

% Eigenvalues
subplot(231)
plot(evals,'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Component number'), ylabel('\lambda')

% Eigenvectors
subplot(232)
plot(evecs(:,1:2),'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Channel'), ylabel('PC weight')
legend({'PC1';'PC2'})

% Original channel modulator
subplot(233)
plot(chanrel,'s-','linew',2,'markerfacecolor','w')
axis square
xlabel('Channel'), ylabel('True channel weights')

% Component time series
subplot(212)
plot(1:N,comp_time_series)
xlabel('Time (a.u.)'), ylabel('Activity')
legend({'PC1';'PC2'})

%% end