%% Single Trial Covariance Distances
%% Load Data
load sampleEEGdata
%% Compute Covariances and Distances
% Compute average of single-trial covariances
% Without storing individual covmats!
covave = zeros( EEG.nbchan );
for triali=1:EEG.trials
    covave = covave + cov( squeeze(EEG.data(:,:,triali))' );
end

% Divide by number of trials
covave = covave / triali;

% Now loop through trials and compute the distance to the average
covdist = zeros(EEG.trials,1);

for triali=1:EEG.trials
    thistrialcov = cov( squeeze(EEG.data(:,:,triali))' );
  
    % Compute Frobenius distance
    covdist(triali) = sqrt( sum(thistrialcov(:) .* covave(:)) );
    % Previous line is the same as this one:
    %covdist(triali) = sqrt( trace(thistrialcov'*covave) );
    
    % Alternative: Euclidean distance (gives similiar results)
    %covdist(triali) = sqrt( sum((thistrialcov(:) - covave(:)).^2) );
end

% Convert to z
covdistz = (covdist-mean(covdist)) / std(covdist);

%% Visual Inspection of Covariance Distances to Average
% Show the covariance distances
figure(1), clf
subplot(2,3,1:2)
plot(covdistz,'ks-','linew',2,'markerfacecolor','w','markersize',12)
xlabel('Trial'), ylabel('Z_{dist}')
title('Z-scored covariance distances')

% Histogram of distances
subplot(233)
histogram(covdistz,10)
xlabel('Distances'), ylabel('Count')
title('Histogram of distances')

%% Pick a Threshold and Reject Trials
% Threshold
thresh = 2.3; % ~.01

% Identify trials that exceed the threshold
toofar = covdistz>thresh;

% Remove those trials from the data
data2 = EEG.data;
data2(:,:,toofar) = [];

%% Show Some Data Before and After Rejection
% Plot time courses
subplot(212), hold on
plot(EEG.times,mean(EEG.data(31,:,:),3),'k','linew',2)
plot(EEG.times,mean(data2(31,:,:),3),'r','linew',2)

% Make the plot look a bit nicer
xlabel('Time (a.u.)')
legend({'Original data';'Trials removed'})
title('Time series before and after covariance cleaning')
zoom on

%% end.