%% Testing normalization
% varying excitatory and suppressive windows to look at the interaction
% between contrast, SOAs, and attention

clear all

%% Simulation 1
% For several (E,S) pairs, vary SOA (50, 100, 200, 400, 800, or more fine
% grained if helpful to smooth the curves), all neutral cues, one contrast.
% You might try a few contrasts to see if the pattern depends on contrast--
% if it does, we'd pick the contrast with strongest effects of (E,S) on SOA.

opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.scaling1 = 2e5;
opt.scaling2 = 3e5;
% opt.aAI = 100;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 400;
opt.tauS1 = 100;
% opt.tR = 2;
% opt.aAV = 2e2;
% opt.AVWeights = [1 0.5];
% opt.distributeVoluntary = 0;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

% loop over model parameters
soas = 1:10;
soa_vals = [100:50:500 800];

windows = [50,400,900];

model_params = combvec(windows,windows,soas);

model_out = nan(2,length(model_params));
for ii=1:length(model_params)
    opt2 = opt;
    rsoa = model_params(3,ii);
    opt2.tauE1 = model_params(1,ii);
    opt2.tauS1 = model_params(2,ii);

    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    model_out(:,ii) = squeeze(perf);
    r1_out(:,:,ii) = p.r1;
end

model_out = reshape(model_out,2,9,10);
r1_out = reshape(r1_out,[size(p.r1) 9 10]);

% plot basic effects
figure
for ii=1:length(windows).^2
    subplot(length(windows),length(windows),ii)
    plot(soa_vals,squeeze(model_out(:,ii,:)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
end

% plot S1 timecourse
figure
for ii=1:length(windows).^2
    subplot(length(windows),length(windows),ii)
    plot(p.tlist,r1_out(:,:,ii))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
end

% sum of r1 layer for each stim
r1_sum = reshape(sum(r1_out([6 12],:,:,:),2),[2 9 10]);

figure
for ii=1:length(windows).^2
    subplot(length(windows),length(windows),ii)
    plot(soa_vals,squeeze(r1_sum(:,ii,:)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
end

%% Simulation 2
% For several (E,S) pairs, at short (250) and long (800) SOAs, vary contrast
% (low / high) independently for T1 and T2, all neutral cues. We want to see
% how (E,S) impacts contrast-dependent suppression, esp at the short SOA.
% I'm thinking of the long SOA as a kind of control where there shouldn't 
% be contrast-dependent suppression, so that would be good to confirm in 
% simulation. If needed we can increase that SOA.

opt = [];
modelClass = [];
rsoa = [4 10]; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64 .64 .32 .32; ...
                     .64 .32 .64 .32];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 400;
opt.tauS1 = 100;
% opt.tR = 2;
% opt.aAV = 2e2;
% opt.AVWeights = [1 0.5];
% opt.distributeVoluntary = 0;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% loop over model parameters
gains = [50,100,400,900];

model_params = combvec(gains,gains);

model_out = nan(2,2,4,length(model_params));
for ii=1:length(model_params)
    opt2 = opt;
    opt2.tauE1 = model_params(1,ii);
    opt2.tauS1 = model_params(2,ii);

    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    model_out(:,:,:,ii) = squeeze(perf);
end

% model_out = model_out-mean(model_out,3);

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:2, squeeze(model_out(1,1,[3 1],ii)))
    plot(1:2, squeeze(model_out(1,1,[4 2],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:2, squeeze(model_out(1,2,[3 1],ii)))
    plot(1:2, squeeze(model_out(1,2,[4 2],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:2, squeeze(model_out(2,1,[2 1],ii)))
    plot(1:2, squeeze(model_out(2,1,[4 3],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:2, squeeze(model_out(2,2,[2 1],ii)))
    plot(1:2, squeeze(model_out(2,2,[4 3],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end

%% Simulation 3
% Same as Simulation 2, except vary attention (attended / neutral / 
% unattended) for T1 and T2 instead of contrast

opt = [];
modelClass = [];
rsoa = [4 10]; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = []; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
% opt.aAI = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 400;
opt.tauS1 = 100;
% opt.tR = 2;
opt.aAV = 2e2;
% opt.AVWeights = [1 0.5];
% opt.distributeVoluntary = 0;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% loop over model parameters
gains = [50,100,400,900];

model_params = combvec(gains,gains);

model_out = nan(2,2,3,length(model_params));
for ii=1:length(model_params)
    opt2 = opt;
    opt2.tauE1 = model_params(1,ii);
    opt2.tauS1 = model_params(2,ii);

    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    model_out(:,:,:,ii) = squeeze(perf);
end

% model_out = model_out-mean(model_out,3);

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:3, squeeze(model_out(1,:,[1 3 2],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 3.2]), xticks(1:3), xticklabels(["A","N","U"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end

figure
for ii=1:length(model_params)
    subplot(length(gains),length(gains),ii), hold on
    plot(1:3, squeeze(model_out(2,:,[2 3 1],ii)))
    title(sprintf('E: %d ms, S: %d ms',model_params(1,ii),model_params(2,ii)))
    xlim([0.8 3.2]), xticks(1:3), xticklabels(["A","N","U"])
    %ylim([0 2.5]), ylabel("Performance (d')")
end
