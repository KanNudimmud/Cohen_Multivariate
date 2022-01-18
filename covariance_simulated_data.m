%% Covariance with Simulated Data
%% Simulate Data with Covariance Structure
% Simulation parameters
N = 1000;     % time points
M =   20;     % channels
nTrials = 50; % number of trials

% Time vector (radian units)
t = linspace(0,6*pi,N);

% Initialize data matrices
data = cell(1,2);

% Relationship across channels (imposing covariance)
chanrel = sin(linspace(0,2*pi,M))';

% Loop over "trials" and generate data
for triali=1:nTrials
    
    % Simulation 1
    data{1}(:,:,triali) = bsxfun(@times,repmat( sin(t),M,1 ),chanrel) + randn(M,N)/1;
    
    % Simulation 2
    data{2}(:,:,triali) = bsxfun(@times,repmat( sin(t+rand*2*pi),M,1 ),chanrel) + randn(M,N)/1;
end

%% Visual Inspection of the Time Series Data
figure(1), clf

for i=1:2
    subplot(1,2,i)
    % Pick a random trial
    trial2plot = ceil(rand*nTrials);
    
    % Show the data
    plot(t,bsxfun(@plus, squeeze(data{i}(:,:,trial2plot)) ,(1:M)'*3 ))
    
    % Code to show trials from one channel
    %plot(t,bsxfun(@plus, squeeze(data{i}(5,:,1:M)) ,(1:M)*3 ))
    
    % Make the plot look a bit nicer
    set(gca,'ytick',[]), axis tight
    xlabel('Time'), ylabel('Channels')
    title([ 'Dataset ' num2str(i) ', trial' num2str(trial2plot) ])
    
    % Code if showing trials instead of channels
    %ylabel('Trials')
    %title([ 'Dataset ' num2str(i) ', channel 5' ])
end

%% Compute Covariance Matrix element-wise
% Initialize covariance matrix
covmat1 = zeros(M);

% Compute the trial average (ERP)
dataave = mean(data{1},3);

% Double-loop over channels and compute dot product scaled by N-1
for chani=1:M
    for chanj=1:M     
        % Mean-center data
        subi = dataave(chani,:) - mean(dataave(chani,:));
        subj = dataave(chanj,:) - mean(dataave(chanj,:));
        
        % Compute covariance
        covmat1(chani,chanj) = sum( subi.*subj ) / (N-1);
    end % end chanj
end % end chani

%% Covariance via Matrix Multiplication
% First mean-center (over time!)
dataaveM = bsxfun(@minus,dataave,mean(dataave,2));

% All pairwise dot products as the matrix times its transpose
covmat2 = dataaveM*dataaveM' / (N-1);

%% Covariance Using MATLAB's cov Function
% Note that cov will mean-center the data
covmat3 = cov(dataave');

% Always check the size to make sure it's the correct orientation:
size(covmat3)

%% Show the Matrices
figure(2), clf
subplot(131), imagesc(covmat1)
axis square, title('Element-wise')
xlabel('Channels'), ylabel('Channels')

subplot(132), imagesc(covmat2)
axis square, title('Matrix multiplication')
xlabel('Channels'), ylabel('Channels')

subplot(133), imagesc(covmat3)
axis square, title('MATLAB cov function')
xlabel('Channels'), ylabel('Channels')

%% Compare the Covariance Matrices for the Two Datasets
[covmat1,covmat2] = deal( zeros(M) );

for triali=1:nTrials
    %% Dataset 1
    % Extract data and mean-center
    tmpdat = squeeze( data{1}(:,:,triali) );
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    % Compute covariance and add to matrix
    covmat1 = covmat1 + tmpdat*tmpdat'/(N-1);
    
    %% Dataset 2
    % Extract data and mean-center
    tmpdat = squeeze( data{2}(:,:,triali) );
    tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
    % Compute covariance and add to matrix
    covmat2 = covmat2 + tmpdat*tmpdat'/(N-1);
end

% Finish the averaging
covmat1 = covmat1 / nTrials;
covmat2 = covmat2 / nTrials;

figure(3), clf
clim = [-1 1]*.7; % colormap
subplot(221), imagesc(covmat1)
axis square, title('Dataset 1, trial covars')
set(gca,'clim',clim,'xtick',[],'ytick',[])

subplot(222), imagesc(covmat2)
axis square, title('Dataset 2, trial covars')
set(gca,'clim',clim,'xtick',[],'ytick',[])

%% Now repeat for average
% Compute average
dataave1 = mean(data{1},3);
dataave2 = mean(data{2},3);

% Covariance for data1
tmpdat  = bsxfun(@minus,dataave1,mean(dataave1,2));
covmat1 = tmpdat*tmpdat' / (N-1);

% Covariance for data2
tmpdat  = bsxfun(@minus,dataave2,mean(dataave2,2));
covmat2 = tmpdat*tmpdat' / (N-1);

subplot(223), imagesc(covmat1)
axis square, title('Dataset 1, average covar')
set(gca,'clim',clim,'xtick',[],'ytick',[])

subplot(224), imagesc(covmat2)
axis square, title('Dataset 2, average covar')
set(gca,'clim',clim,'xtick',[],'ytick',[])
set(gca,'clim',[-1 1]*.05)

%% end