function testTemporalWindows

opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

% H/H, H/L, L/H, L/L
% opt.stimContrasts = [.64 .64 .16 .16; ...
%                      .64 .16 .64 .16];
opt.stimContrasts = [.64; .64];
% opt.scaling1 = 3e4;
% opt.scaling2 = 3e4;
% opt.aAI = 0;
% opt.tR = 2;
% opt.distributeVoluntary = 0;
% opt.AVWeights = [1 0.8];

% opt.tau1 = 2;
% opt.tau2 = 2;

opt.display.plotTS = 1; % plot the time series for each simulation

% loop over excitatory and inhibitory window params
tauE_list = 0:50:1000;
tauS_list = 0:50:1000;

tau_list = CombVec(tauE_list,tauS_list);
perf_out = nan(length(tau_list),2);

parfor ii=1:length(tau_list)
    opt2 = opt;
    opt2.tauE1 = tau_list(1,ii);
    opt2.tauS1 = tau_list(2,ii);
    [~,~,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    perf_out(ii,:) = perf;
end

perf_out = reshape(perf_out,[length(tauE_list),length(tauS_list),2]);

figure
imagesc(tauE_list,tauS_list,perf_out(:,:,1))
colorbar
axis square

figure
imagesc(tauE_list,tauS_list,perf_out(:,:,2))
colorbar
axis square
