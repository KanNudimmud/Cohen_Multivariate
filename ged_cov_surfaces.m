%% Source separation with GED
% GED intuition with covariance surfaces
%% Quadratic Forms of S, R, and S/R
% Create two little matrices
S = [1 2; 2 9];
R = [2 -1; -1 4];

% Weights along each dimension
wi = -2:.1:2;

% Initialize
[quadformS,quadformR] = deal( zeros(length(wi)) );
quadformG = nan(length(wi));

% Fudge-factor for preventing the surface from running off
fudge = .1;

for i=1:length(wi)
    for j=1:length(wi)
        % Define the weights for this iteration
        w  = [wi(i) wi(j)]';
        nf = w'*w;
        
        % Quadratic form of matrix S
        quadformS(i,j) = w'*S*w / nf;
        
        % QF of matrix R
        quadformR(i,j) = w'*R*w / nf;
        
        % QF of ratio
        quadformG(i,j) = (w'*S*w+fudge) / (w'*R*w+fudge);
    end
end

% Eigenvectors of S and R
[evecsS,evals] = eig(S);
[evecsR,evals] = eig(R);

% Eigenvectors from GED
[evecsGED,evals] = eig(S,R);

%% Visualize the Eigenvectors
figure(1), clf
subplot(211), hold on

plot([0 evecsS(1,1)],[0 evecsS(2,1)],'r','linew',3)
plot([0 evecsS(1,2)],[0 evecsS(2,2)],'r','linew',3)

plot([0 evecsR(1,1)],[0 evecsR(2,1)],'m','linew',3)
plot([0 evecsR(1,2)],[0 evecsR(2,2)],'m','linew',3)

plot([0 evecsGED(1,1)],[0 evecsGED(2,1)],'k','linew',3)
plot([0 evecsGED(1,2)],[0 evecsGED(2,2)],'k','linew',3)

axis([-1 1 -1 1]*1)
grid on, axis square
legend({'S_1','S_2','R_1','R_2','GED_1','GED_2'})
title('Eigenvectors')

subplot(234)
imagesc( evecsS'*evecsS )
axis square
set(gca,'clim',[0 1])
title('W_S^TW_S')

subplot(235)
imagesc( evecsR'*evecsR )
axis square
set(gca,'clim',[0 1])
title('W_R^TW_R')

subplot(236)
imagesc( evecsGED'*evecsGED )
axis square
set(gca,'clim',[0 1])
title('W_G^TW_G')

%% Some Nice Butterflies
figure(2), clf

subplot(231), imagesc(S), axis square, title('Matrix S')
set(gca,'xtick',[],'ytick',[])

subplot(232), imagesc(R), axis square, title('Matrix R')
set(gca,'xtick',[],'ytick',[])

subplot(233), imagesc(inv(R)*S), axis square, title('Matrix R^{-1}S')
set(gca,'xtick',[],'ytick',[])

subplot(234)
surf(wi,wi,quadformS'), shading interp
title('Quadratic form of S')
xlabel('W1'), ylabel('W2'), zlabel('Energy'), axis square

subplot(235)
surf(wi,wi,quadformR'), shading interp
title('Quadratic form of R')
xlabel('W1'), ylabel('W2'), zlabel('Energy'), axis square

subplot(236)
surf(wi,wi,quadformG'), shading interp
title('Quadratic form of S/R')
xlabel('W1'), ylabel('W2'), zlabel('Energy ratio'), axis square

%% end.