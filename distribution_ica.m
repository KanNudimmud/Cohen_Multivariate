%% Independent components analysis
% Data distributions and ICA
%% 
% A number
N = 1000;

% A non-Gaussian distribution
dist1 = rand(N,1);

% Another non-Gaussian distribution
dist2 = rand(N,1).^2;

% Their separate and combined distributions
figure(1), clf
subplot(221)
histogram(dist1,100)
title('Distribution 1')

subplot(222)
histogram(dist2,100)
title('Distribution 1')

subplot(212)
histogram(dist1+dist2,100)
title('Distribution 1+2')

%% ICA on a trivial example (unmixed data)
clear data

% Wwo non-Gaussian distributions
data(1,:) = rand(1,N);
data(2,:) = rand(1,N).^4;

% ICA and scores
b = jader(data);
iscores = b*data;

% Plot distributions
figure(2), clf

% Data 1
subplot(221), histogram(data(1,:),100)
title('Data 1')

% Data 2
subplot(222), histogram(data(2,:),100)
title('Data 2')

% IC 1
subplot(223), histogram(iscores(1,:),100)
title('IC 1')

% IC 2
subplot(224), histogram(iscores(2,:),100)
title('IC 2')

% Plot data as a function of ICs
figure(3), clf
subplot(121)
plot(data(1,:),iscores(1,:),'o')
xlabel('Data'), ylabel('IC1 scores')
axis square

subplot(122)
plot(data(2,:),iscores(2,:),'o')
xlabel('Data'), ylabel('IC2 scores')
axis square

%% Try again with mixed data
clear data

% Two non-Gaussian distributions
dataO(1,:) = rand(1,N);
dataO(2,:) = rand(1,N).^4;

% Then mix them into a new dataset
data(1,:) = .3*dataO(1,:) + .7*dataO(2,:);
data(2,:) = .7*dataO(1,:) + .3*dataO(2,:);

% ICA and scores (using the mixed data!)
b = jader(data);
iscores = b*data;

% Plot distributions
figure(4), clf

% Original data 1
subplot(321), histogram(dataO(1,:),100)
title('O.D. 1')

% Original data 2
subplot(322), histogram(dataO(2,:),100)
title('O.D. 2')

% Mixed data 1
subplot(323), histogram(data(1,:),100)
title('Mixed data 1')

% Mixed data 2
subplot(324), histogram(data(2,:),100)
title('Mixed data 2')

% IC 1
subplot(325), histogram(iscores(1,:),100)
title('IC 1')

% IC 2
subplot(326), histogram(iscores(2,:),100)
title('IC 2')

% Plot data as a function of ICs
figure(5), clf
subplot(121)
plot(data(1,:),iscores(1,:),'o')
xlabel('Data'), ylabel('IC1 scores')
axis square

subplot(122)
plot(data(2,:),iscores(2,:),'o')
xlabel('Data'), ylabel('IC2 scores')
axis square

%% Distribution shapes of non-stationary sine waves
%%% The purpose of this cell is for you to see the signal value distributions 
%   of sine waves with varying degrees of non-stationarities.
%   Note: The code is hard-coded to five sine waves.
% Some parameters
pnts  = 6000;
srate = 1000;

% Frequencies
hz = linspace(0,srate,pnts);

% Specify a range of FWHM
minFWHM = .1; % in Hz
maxFWHM = 10; % also in Hz.

% Specify range of FWHM
fwhm = logspace(log10(.1),log10(10),5);
peakfreq = 8;

figure(6), clf

for fwhmi=1:5
    % Create frequency-domain Gaussian
    s  = fwhm(fwhmi)*(2*pi-1)/(4*pi); % normalized width
    x  = hz-peakfreq;          % shifted frequencies
    fg = exp(-.5*(x/s).^2);    % gaussian
    
    % Random Fourier coefficients and taper with the Gaussian
    fc = rand(1,pnts) .* exp(1i*2*pi*rand(1,pnts));
    fc = fc .* fg;
    
    % Create the time-domain signal, and normalize
    sig = real( ifft(fc) );
    sig = (sig-mean(sig)) / std(sig);
    
    %%% Plot the power spectrum
    subplot(5,3,fwhmi*3-2)
    plot(hz,abs(fc).^2)
    set(gca,'xlim',[0 peakfreq*4],'xtick',[],'ytick',[]), box off
    title([ 'FWHM = ' num2str(round(fwhm(fwhmi),2)) ])
    if fwhmi==5, xlabel('Frequency (Hz)'), end
    
    %%% Plot the time-domain signal
    subplot(5,3,fwhmi*3-1)
    plot(sig), box off
    set(gca,'xtick',[],'ytick',[])
    if fwhmi==5, xlabel('Time'), end
    
    %%% Plot the histogram
    subplot(5,3,fwhmi*3)
    [yh,xh]=hist(sig,90);
    plot(xh,yh)
    set(gca,'xlim',[-1 1]*5,'xtick',[],'ytick',[]), box off
    if fwhmi==5, xlabel('Data value'), end
end

%% end.