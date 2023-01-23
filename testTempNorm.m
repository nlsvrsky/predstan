
opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.tau1_t = 10;
opt.tau2_t = 60;

opt.display.plotTS = 1; % plot the time series for each simulation

runModel(opt, modelClass, rsoa, rseq, rcond);
