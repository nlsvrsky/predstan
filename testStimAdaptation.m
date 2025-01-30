%% test stimulus adaptation

opt = [];
modelClass = [];
rsoa = 500; % SOA = 250 ms (see runModel)
rseq = 3; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

opt.stimContrasts = [.64; .64];
opt.stimDur = 30;
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 100;
opt.tauS1 = 50;
% opt.tau1 = 2;
% opt.sigma1 = .01;

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

% opt.tauD = 500;
% opt.tR = 2;
% opt.aAV = 1e2;
% opt.AVWeights = [1 1];
% opt.AVNeutralT1Weight = 0.25;
% opt.distributeVoluntary = 1;

opt.display.plotTS = 1; % plot the time series for each simulation
opt.display.plotPerf = 0;

[~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
