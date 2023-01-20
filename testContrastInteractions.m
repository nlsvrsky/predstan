% model testing
opt = [];
% H/H, H/L, L/H, L/L
opt.stimContrasts = [.64 .64 .16 .16; ...
                     .64 .16 .64 .16];

modelClass = [];
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
