%% find gain parameters that maximize performance for T1 and T2
p = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = []; % cueT1, cueT2

% p.tauE1 = 20;
% p.tauS1 = 50;
% p.tau2 = 2;
p.sigma1 = .1;

% p.tauE2 = 50;
% p.tauS2 = 100;
% p.tau2 = 2;
p.sigma2 = .1;

% p.tauEAV = 50;
% p.tauSAV = 50;
% p.tauAV = 2;
% p.sigmaA = 100;

% p.tR = 1e5;
% p.distributeVoluntary = 0;
p.aAV = 100;
% p.AVWeights = [1 1];

% p.AVOnset = 10;
% p.AVDur = 40;

p.display.plotTS = 0; % plot the time series for each simulation

% loop over contrasts
contrMin = 0.01;
contrMax = 0.95;
contrN = 10;
fullContrast = exp(linspace(log(contrMin),log(contrMax),contrN));

contrPerf = nan(2,3,length(fullContrast));
sortedPerf = nan(size(contrPerf));

parfor cc=1:length(fullContrast)
    opt = p;
    opt.stimContrasts = [fullContrast(cc); fullContrast(cc)];

    [~,~,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    contrPerf(:,:,cc) = squeeze(perf);
end

sortedPerf(1,1,:) = contrPerf(1,1,:); % T1 valid
sortedPerf(1,2,:) = contrPerf(1,3,:); % T1 neutral
sortedPerf(1,3,:) = contrPerf(1,2,:); % T1 invalid

sortedPerf(2,1,:) = contrPerf(2,2,:); % T2 valid
sortedPerf(2,2,:) = contrPerf(2,3,:); % T2 neutral
sortedPerf(2,3,:) = contrPerf(2,1,:); % T2 invalid

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
