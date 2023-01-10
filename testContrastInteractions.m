% model testing
opt = [];
% H/H, H/L, L/H, L/L
opt.stimContrasts = [.64 .64 .32 .32; ...
                     .64 .32 .64 .32];

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

conds = ["high/high","high/low","low/high","low/low"];
figure
for ii=1:4
    subplot(1,4,ii)
    plot(1:3, sortedPerf(:,:,ii))
    title(sprintf("Contrasts: %s",conds(ii)))
    xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
    ylabel("Performance (d')"), ylim([0 3])
end
legend(["T1","T2"])
