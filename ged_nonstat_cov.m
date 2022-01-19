%% Source separation with GED
% GED with nonstationary covariances
%% Simulate Data with Covariance Structure
% Simulation parameters
N = 500; % time points per covariance matrix
M =  20; % channels

% Time vector (radian units)
t = linspace(0,6*pi,N);

% Relationship across channels (imposing covariance)
chanrel1 = randn(M,1);
chanrel2 = sin(linspace(0,2*pi,M))';
chanrel3 = exp( -linspace(-2,2,M).^2 )';

% All parts of data
data1 = bsxfun(@times,repmat( sin(t+rand*2*pi),M,1 ),chanrel1) + randn(M,N);
data2 = bsxfun(@times,repmat( sin(t+rand*2*pi),M,1 ),chanrel2) + randn(M,N);
data3 = bsxfun(@times,repmat( sin(t+rand*2*pi),M,1 ),chanrel3) + randn(M,N);

% Show the data
figure(1), clf, hold on
plot(bsxfun(@plus,[data1 data2 data3]',(1:M)*3))
axis tight
plot(repmat([N 2*N],2,1),repmat(get(gca,'ylim')',1,2),'w--','linew',3)
set(gca,'ytick',[])
xlabel('Time') ,ylabel('Channels')

%% Covariance Matrices of the Three Windows
% Covariances in a loop with eval
for i=1:3
    eval([ 'tmpdat = bsxfun(@minus,data' num2str(i) ',mean(data' num2str(i) ',2)); ' ]);
    eval([ 'covm' num2str(i) '  = tmpdat*tmpdat'' / (N-1);' ]);
end

% Covariance of the 2nd and 3rd data chunks
tmpdat = [data2 data3];
tmpdat = bsxfun(@minus,tmpdat,mean(tmpdat,2));
covm23 = tmpdat*tmpdat'/(2*N-1);

% Plot the covariance matrices
clim = [-1 1]*.5;
figure(2), clf
for i=1:4
    subplot(2,2,i)
    
    if i<4
        eval([ 'imagesc(covm' num2str(i) ')' ])
        title([ 'Data ' num2str(i) ])
    else
        imagesc(covm23)
        title('Data 2&3')
    end
    axis square
    set(gca,'clim',clim,'xtick',[],'ytick',[])
    xlabel('Channels'), ylabel('Channels')
end

%% GED Relative to the First Time Window
[V2,D2]   = eig( covm2,covm1 );
[D2,sidx] = sort( diag(D2),'descend' );
V2 = V2(:,sidx);

[V3,D3]   = eig( covm3,covm1 );
[D3,sidx] = sort( diag(D3),'descend' );
V3 = V3(:,sidx);

[V23,D23]  = eig( covm23,covm1 );
[D23,sidx] = sort( diag(D23),'descend' );
V23 = V23(:,sidx);

%% Plot the Filter Forward Models
%%% Axis limits to make plots comparable
% maximum of all eigenvalues for plotting
maxy = max([ D2;D3;D23 ])*1.1;
ylim = [-1 1]*2;

figure(3), clf

%%% Period 2 vs. baseline
subplot(231)
plot(D2,'ks-','markerfacecolor','w','markersize',10)
xlabel('Component number'), ylabel('\lambda')
set(gca,'ylim',[0 maxy])
axis square
title('GED: 2 vs 1')

subplot(234), hold on
plot(V2(:,1)'*covm2,'ks-','markerfacecolor','w','markersize',10)
plot(chanrel2,'ro-','markerfacecolor','w','markersize',10)
axis square, title('Filter forward model')
set(gca,'ylim',ylim)
xlabel('Channel'), ylabel('Weight')
legend({'Filter';'Simulated'})

%%% Period 3 vs. baseline
subplot(232)
plot(D3,'ks-','markerfacecolor','w','markersize',10)
xlabel('Component number'), ylabel('\lambda')
set(gca,'ylim',[0 maxy])
axis square
title('GED: 3 vs. 1')

subplot(235), hold on
plot(V3(:,1)'*covm3,'ks-','markerfacecolor','w','markersize',10)
plot(chanrel3,'ro-','markerfacecolor','w','markersize',10)
set(gca,'ylim',ylim)
xlabel('Channel'), ylabel('Weight')
axis square, title('Filter forward model')

%%% Periods 2&3 vs. baseline
subplot(233)
plot(D23,'ks-','markerfacecolor','w','markersize',10)
xlabel('Component number'), ylabel('\lambda')
set(gca,'ylim',[0 maxy])
axis square
title('GED: 2&3 vs 1')

subplot(236), hold on
plot(V23(:,1)'*covm23,'ks-','markerfacecolor','w','markersize',10)
plot(V23(:,2)'*covm23,'bs-','markerfacecolor','w','markersize',10)
xlabel('Channel'), ylabel('Weight')
set(gca,'ylim',ylim)
legend({'Filter 1';'Filter 2'})
axis square, title('Filter forward model')

%% end.