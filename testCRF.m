% get contrast response function from DNMA
%clear all; close all;

contrMin = 0.01;
contrMax = 0.64;
contrN = 20;
fullContrast = exp(linspace(log(contrMin),log(contrMax),contrN));

modelClass = [];
opt = [];
opt.stimContrasts = [fullContrast; fullContrast];
opt.sigma1 = 0.1;
opt.sigma2 = 0.01;
opt.sigmaD = 0.7;

opt.scaling1 = 4e4;
opt.scaling2 = 4e4;

rsoa = 4; % SOA = 250 ms
rseq = []; % default orientation sequence
rcond = []; % cueT1, cueT2, cueN

opt.display.plotTS = 0; % plot the time series for each simulation

[~,~,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);

perf = squeeze(perf);
sortedPerf = nan(size(perf));

sortedPerf(1,1,:) = perf(1,1,:); % T1 valid
sortedPerf(1,2,:) = perf(1,3,:); % T1 neutral
sortedPerf(1,3,:) = perf(1,2,:); % T1 invalid

sortedPerf(2,1,:) = perf(2,2,:); % T2 valid
sortedPerf(2,2,:) = perf(2,3,:); % T2 neutral
sortedPerf(2,3,:) = perf(2,1,:); % T2 invalid

sortedPerf = normcdf(sortedPerf./sqrt(2));

% plot CRFs
figure
subplot(1,2,1)
semilogx(fullContrast, squeeze(sortedPerf(1,:,:)))
title("T1")
xlabel("Stimulus contrast"), %xlim([0.1 1]), xticks([.1 .2:.2:1])
ylabel("Performance (d')"), %ylim([0 3])

subplot(1,2,2)
semilogx(fullContrast, squeeze(sortedPerf(2,:,:)))
title("T2")
xlabel("Stimulus contrast"), %xlim([0.1 1]), xticks([.1 .2:.2:1])
ylabel("Performance (d')"), %ylim([0 3])
legend(["Valid","Neutral","Invalid"])
