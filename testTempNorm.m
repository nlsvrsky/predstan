
opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

% H/H, H/L, L/H, L/L
% opt.stimContrasts = [.9; .9];
% opt.stimContrasts = [.64 .64 .16 .16; ...
%                      .64 .16 .64 .16];
% opt.stimContrasts = [.64 .64 .32 .32; ...
%                      .64 .32 .64 .32];
opt.stimContrasts = [.64; .64];
% opt.scaling1 = 4e5;
% opt.scaling2 = 5e5;
opt.aAI = 0;

opt.tauE1 = 800;
opt.tauS1 = 800;
% opt.tau1 = 2;
% opt.sigma1 = 8;

opt.tauE2 = 0;
opt.tauS2 = 0;
% opt.tau2 = 5000;
% opt.sigma2 = 10;

% opt.tauEAV = 50;
% opt.tauSAV = 50;
% opt.tauAV = 2;
% opt.sigmaA = 100;

% opt.AVOnset = -20;
% opt.AVDur = 60;

% opt.tauD = 500;
% opt.tR = 1e5;
% opt.aAV = 1000;
% opt.AVWeights = [.8 1];
% opt.AVNeutralT1Weight = 0.25;
% opt.distributeVoluntary = 0;

opt.display.plotTS = 1; % plot the time series for each simulation

[~,~,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);

if isempty(rcond)
    perf = squeeze(perf);
    sortedPerf = nan(size(perf));
    
    sortedPerf(1,1,:) = perf(1,1,:); % T1 valid
    sortedPerf(1,2,:) = perf(1,3,:); % T1 neutral
    sortedPerf(1,3,:) = perf(1,2,:); % T1 invalid
    
    sortedPerf(2,1,:) = perf(2,2,:); % T2 valid
    sortedPerf(2,2,:) = perf(2,3,:); % T2 neutral
    sortedPerf(2,3,:) = perf(2,1,:); % T2 invalid
    
    conds = ["T1 high/T2 high","T1 high/T2 low","T1 low/T2 high","T1 low/T2 low"];
    figure
    for ii=1:4
        subplot(1,4,ii)
        plot(1:3, sortedPerf(:,:,ii))
        title(sprintf("Contrasts: %s",conds(ii)))
        xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
        ylabel("Performance (d')"), ylim([0 3])
    end
    legend(["T1","T2"])
    
    % neutral condition only
    figure
    subplot(1,2,1), hold on, title("T1")
    plot(1:2, squeeze(sortedPerf(1,2,[3 1]))) % T1 with T2 high contrast
    plot(1:2, squeeze(sortedPerf(1,2,[4 2]))) % T1 with T2 low contrast
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    ylim([0 3]), ylabel("Performance (d')")
    legend(["T2 high","T2 low"])
    
    subplot(1,2,2), hold on, title("T2")
    plot(1:2, squeeze(sortedPerf(2,2,[2 1]))) % T2 with T1 high contrast
    plot(1:2, squeeze(sortedPerf(2,2,[4 3]))) % T2 with T1 low contrast
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
    ylim([0 3]), ylabel("Performance (d')")
    legend(["T1 high","T1 low"])

elseif rcond==3
    figure
    subplot(1,2,1), hold on, title("T1")
    plot(1:2, perf(1,[3 1]))
    plot(1:2, perf(1,[4 2]))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    ylim([0 2.5]), ylabel("Performance (d')")
    legend(["T2 high","T2 low"])
    
    subplot(1,2,2), hold on, title("T2")
    plot(1:2, perf(2,[2 1])) % T2 with T1 high contrast
    plot(1:2, perf(2,[4 3])) % T2 with T1 low contrast
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
    ylim([0 2.5]), ylabel("Performance (d')")
    legend(["T1 high","T1 low"])
end
