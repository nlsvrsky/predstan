% get contrast response function from DNMA
%clear all; close all;

contrMin = 0.01;
contrMax = 0.95;
contrN = 10;
fullContrast = exp(linspace(log(contrMin),log(contrMax),contrN));

modelClass = [];
opt = [];
% opt.stimContrasts = [fullContrast; fullContrast];
opt.sigma1 = 5;
opt.sigma2 = 0.1;
opt.sigmaD = 0.7;

opt.scaling1 = 1e5;
opt.scaling2 = 1e5;

opt.aAV = 250;

opt.AVOnset = -150;
opt.AVDur = 250;

rsoa = 7; % SOA = 250 ms
rseq = []; % default orientation sequence
rcond = []; % cueT1, cueT2, cueN

opt.display.plotTS = 0; % plot the time series for each simulation

perfs = nan([2 3 contrN]);
sortedPerf = nan(size(perfs));

parfor ii=1:contrN
    thisOpt = opt;
    thisOpt.stimContrasts = [fullContrast(ii); fullContrast(ii)];
    [~,~,perf] = runModel(thisOpt, modelClass, rsoa, rseq, rcond);
    perfs(:,:,ii) = squeeze(perf);
end

sortedPerf(1,1,:) = perfs(1,1,:); % T1 valid
sortedPerf(1,2,:) = perfs(1,3,:); % T1 neutral
sortedPerf(1,3,:) = perfs(1,2,:); % T1 invalid

sortedPerf(2,1,:) = perfs(2,2,:); % T2 valid
sortedPerf(2,2,:) = perfs(2,3,:); % T2 neutral
sortedPerf(2,3,:) = perfs(2,1,:); % T2 invalid

% sortedPerf = normcdf(sortedPerf./sqrt(2));

% plot CRFs
figure
subplot(1,2,1)
semilogx(fullContrast, squeeze(sortedPerf(1,:,:)))
title("T1")
xlabel("Stimulus contrast"), %xlim([0.1 1]), xticks([.1 .2:.2:1])
ylabel("Performance (d')")%, ylim([0 3])

subplot(1,2,2)
semilogx(fullContrast, squeeze(sortedPerf(2,:,:)))
title("T2")
xlabel("Stimulus contrast"), %xlim([0.1 1]), xticks([.1 .2:.2:1])
ylabel("Performance (d')")%, ylim([0 3])
legend(["Valid","Neutral","Invalid"])
