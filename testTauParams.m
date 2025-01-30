
opt = [];
modelClass = [];
rsoa = 250; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

opt.scaling1 = 6e5;
opt.scaling2 = 6e5;
opt.aAI = 0;
opt.aAV = 0;

opt.eScale = 1;
opt.sScale = 1;

opt.tauE1 = 0;
opt.tauS1 = 400;
% opt.tau1 = 2;
opt.sigma1 = .1;

% opt.tauE2 = 0;
% opt.tauS2 = 0;
% opt.tau2 = 5000;
% opt.sigma2 = 10;

% opt.tauEAV = 50;
% opt.tauSAV = 50;
% opt.tauAV = 2;
% opt.sigmaA = 100;

% opt.AVOnset = 0;
% opt.AVDur = 60;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

tau_list = 0:50:1000;
contr_list = [.64 .64 .32; .64 .32 .64];

params = combvec(tau_list,tau_list,contr_list);

perf_out = nan(2,length(params));
parfor ii=1:length(params)
    % fprintf('%d/%d\n',ii,length(params));
    opt2 = opt;
    opt2.stimContrasts = params(3:4,ii);
    opt2.tauE1 = params(1,ii);
    opt2.tauS1 = params(2,ii);
    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    perf_out(:,ii) = perf;
end

perf_out = reshape(perf_out,[2 length(tau_list) length(tau_list) length(contr_list)]);

perf_T1 = squeeze((perf_out(1,:,:,2)-perf_out(1,:,:,1))./(perf_out(1,:,:,2)+perf_out(1,:,:,1)));
perf_T2 = squeeze((perf_out(2,:,:,3)-perf_out(2,:,:,1))./(perf_out(2,:,:,3)+perf_out(2,:,:,1)));

figure
subplot(131), imagesc(tau_list,tau_list,perf_T1'), axis square, clim([0 max(perf_T1,[],'all')])
subplot(132), imagesc(tau_list,tau_list,perf_T2'), axis square, clim([0 max(perf_T2,[],'all')])
subplot(133), imagesc(tau_list,tau_list,perf_T1'.*perf_T2'), axis square, clim([0 max(perf_T1.*perf_T2,[],'all')])
