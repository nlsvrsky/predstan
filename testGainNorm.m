%% Testing gain parameters for model with normalization
%  varying gain applied to each stimulus to determine the optimal
%  combination for different cueing conditions

clear all

opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 1; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.scaling1 = 6e5;
opt.scaling2 = 7e5;
opt.aAI = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 400;
opt.tauS1 = 100;
opt.tR = 2;
opt.aAV = 1e3;
opt.AVWeights = [1 0.5];
opt.distributeVoluntary = 0;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% loop over model parameters
soas = 1:10;
soa_vals = [100:50:500 800];

weights = 0:0.1:1;

model_params = combvec(weights,weights,soas);

model_out = nan(2,length(model_params));
parfor ii=1:length(model_params)
    opt2 = opt;
    rsoa = model_params(3,ii);
    opt2.AVWeights = model_params(1:2,ii);

    [~,~,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    model_out(:,ii) = squeeze(perf);
end

model_out = reshape(model_out,[2 length(weights) length(weights) length(soas)]);

figure(1)
for ii=1:2, for jj=1:length(soas)
    subplot(length(soas),2,(jj-1)*2+ii)
    imagesc(weights,weights,squeeze(model_out(ii,:,:,jj)))
    title(sprintf("T%d, %d ms",ii,soa_vals(jj)))
    clim([min(model_out(ii,:,:,:),[],'all') max(model_out(ii,:,:,:),[],'all')])
end, end

figure(2)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.5;.5].*model_out(:,:,:,jj))))
    title(sprintf("Neutral, %d ms",soa_vals(jj)))
%     clim([min(sum([.5;.5].*model_out),[],'all') max(sum([.5;.5].*model_out),[],'all')])
end

figure(3)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.75;.25].*model_out(:,:,:,jj))))
    title(sprintf("Cue T1, %d ms",soa_vals(jj)))
%     clim([min(sum([.75;.25].*model_out),[],'all') max(sum([.75;.25].*model_out),[],'all')])
end

figure(4)
for jj=1:length(soas)
    subplot(2,5,jj)
    imagesc(weights,weights,squeeze(sum([.25;.75].*model_out(:,:,:,jj))))
    title(sprintf("Cue T2, %d ms",soa_vals(jj)))
%     clim([min(sum([.25;.75].*model_out),[],'all') max(sum([.25;.75].*model_out),[],'all')])
end

% T1 fixed at max, function of SOA and T2 gain
figure(5), subplot(121)
plot(weights,squeeze(sum([.75;.25].*model_out(:,end,:,:))))
title('Cue T1, with T1 at max gain')
xlabel('T2 Weight'), ylabel("Expected d'")

% T2 fixed at max, function of SOA and T1 gain
figure(5), subplot(122)
plot(weights,squeeze(sum([.25;.75].*model_out(:,:,end,:))))
title('Cue T2, with T2 at max gain')
xlabel('T1 Weight'), ylabel("Expected d'")
