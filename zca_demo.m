%% Source separation with GED
% ZCA demo
%% 
d = linspace(.2,2,14);

figure(1), clf
plot(d,d.^(-1/2),'ks-','linew',2,'markerfacecolor','w','markersize',13)

hold on
plot([d(1) 1],[1 1],'k--')
yy = get(gca,'ylim');
plot([1 1],[yy(1) 1],'k--')

axis square
xlabel('d')
ylabel('d^{-1/2}')

%% Generate 2D data
% Number of data points
N = 1000;

% Part 1 of the data
x1 = [ 1*randn(N/2,1) .1*randn(N/2,1) ];

% Rotation matrix
th = pi/4;
R1 = [ cos(th) -sin(th) ;
       sin(th)  cos(th) ];

% Part 2 of the data
x2 = [ 1*randn(N/2,1) .1*randn(N/2,1) ];

% Note the different rotation angle
th = pi;
R2 = [ cos(th) -sin(th) ;
       sin(th)  cos(th) ];

% Combine and rotate data
y = [x1*R1; x2*R2];

%% Eigendecomposition for PCA and ZCA
% Eig of y
y = bsxfun(@minus,y,mean(y,1));
covmat = y'*y / (length(y)-1);
[V,D] = eig(covmat);

% Component data
c = y*V;

% ZCA
yz = ( V*D^(-1/2)*V'*y' )';

% PCA of ZCA data
[Vz,Dz] = eig( yz'*yz );
cz = yz*Vz;

%% Plotting
figure(1), clf

% Original data
subplot(221), hold on
plot(y(:,1),y(:,2),'ko','markerfacecolor','m')
plot([0 V(1,1)],[0 V(2,1)],'k','linew',3)
plot([0 V(1,2)],[0 V(2,2)],'k','linew',3)
axis square
axis([-1 1 -1 1]*max(abs(y(:))))
xlabel('y_1'), ylabel('y_2')
title('Data')

% Component projections
subplot(222)
plot(c(:,1),c(:,2),'ko','markerfacecolor','m')
axis square
axis([-1 1 -1 1]*max(abs(y(:))))
xlabel('pc_1'), ylabel('pc_2')
title('PCA of data')

% Whitened data
subplot(223), hold on
plot(yz(:,1),yz(:,2),'ko','markerfacecolor','m')
plot(yz(1,1),yz(1,2),'rs','markerfacecolor','y','markersize',15)
axis square
axis([-1 1 -1 1]*max(abs(yz(:))))
xlabel('y_{z1}'), ylabel('y_{z2}')
title('Whitened')
plot([0 Vz(1,1)],[0 Vz(2,1)],'k','linew',3)
plot([0 Vz(1,2)],[0 Vz(2,2)],'k','linew',3)

subplot(224), hold on
plot(cz(:,1),cz(:,2),'ko','markerfacecolor','m')
plot(cz(1,1),cz(1,2),'rs','markerfacecolor','y','markersize',15)
axis square
axis([-1 1 -1 1]*max(abs(yz(:))))
xlabel('pc_1'), ylabel('pc_2')
title('PCA of whitened')

%% end.