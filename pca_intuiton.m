%% PCA intuition with scatter plots and covariance surfaces
%% PCA on simulated data
% Data
x = [ 1*randn(1000,1) .4*randn(1000,1) ];

% Rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th); sin(th) cos(th) ];

% Rotate data to induce correlation
y = x*R1;

% Plot the data in a scatter plot
figure(1), clf
subplot(131)
plot(y(:,1),y(:,2),'m.','markersize',17)

% Make the plot look nicer
set(gca,'xlim',[-1 1]*max(y(:)),'ylim',[-1 1]*max(y(:)))
xlabel('x-axis'), ylabel('y-axis')
axis square

%% PCA via eigendecomposition
% Mean-center
y = y - mean(y,1);

% Covariance matrix
covmat = y'*y / (length(y)-1);

% Eigendecomposition
[evecs,evals] = eig(covmat);

% Plot the eigenvectors on the data
hold on
plot([0 evecs(1,1)]*evals(1),[0 evecs(2,1)]*evals(1),'k','linew',3)
plot([0 evecs(1,2)]*evals(end),[0 evecs(2,2)]*evals(end),'k','linew',3)

% Show the covariance matrix
subplot(132)
imagesc(covmat), axis square
set(gca,'clim',[-1 1],'xtick',1:2,'ytick',1:2)
xlabel('Channels'), ylabel('Channels')

%% Compute quadratic form
% Weights along each dimension
xi = -2:.1:2;

quadform = zeros(length(xi));
for i=1:length(xi)
    for j=1:length(xi)
        % Define vector
        x = [xi(i) xi(j)]';
        
        % QF
        quadform(i,j) = x'*covmat*x / (x'*x);
    end
end

% Fill in missing point with 0
quadform(~isfinite(quadform)) = 0;

%% Visualize the quadratic form surface
subplot(133)
surf(xi,xi,quadform'), shading interp
xlabel('W1'), ylabel('W2'), zlabel('energy')
rotate3d on, axis square

% Eigenvectors
hold on
plot3([0 evecs(1,1)]*2,[0 evecs(2,1)]*2,[1 1]*max(quadform(:)),'w','linew',3)
plot3([0 evecs(1,2)]*2,[0 evecs(2,2)]*2,[1 1]*max(quadform(:)),'m','linew',3)

%% end.