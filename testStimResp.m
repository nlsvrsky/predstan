% test stimulus responses

opt = [];
modelClass = [];
rsoa = 10; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

% H/H, H/L, L/H, L/L
% opt.stimContrasts = [.9; .9];
% opt.stimContrasts = [.64 .64 .16 .16; ...
%                      .64 .16 .64 .16];
% opt.stimContrasts = [.64 .64 .32 .32; ...
%                      .64 .32 .64 .32];
opt.stimContrasts = [.64; 0];
opt.stimDur = 500;
opt.scaling1 = 5e5;
opt.scaling2 = 6e5;
% opt.aAI = 1e2;

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
opt.aAV = 2e2;
% opt.AVWeights = [1 1];
% opt.AVNeutralT1Weight = 0.25;
% opt.distributeVoluntary = 0;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

tauE_list = [50 100:100:800];

for ii=1:length(tauE_list)
    opt.tauE1 = tauE_list(ii);
    opt.tauS1 = 0;
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_E(:,:,ii) = p.r1;

    opt.tauS1 = tauE_list(ii);
    opt.tauE1 = 0;
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_S(:,:,ii) = p.r1;
    d1_S(:,:,ii) = p.d1;
    s1_S(:,:,ii) = p.s1;
    f1_S(:,:,ii) = p.f1;

    rai(:,:,ii) = p.rai;
    dai(:,:,ii) = p.dai;
    sai(:,:,ii) = p.sai;
    fai(:,:,ii) = p.fai;
end

figure
subplot(121)
plot(p.tlist,squeeze(r1_E(6,:,:)))

subplot(122)
plot(p.tlist,squeeze(r1_S(6,:,:)))

%% subadditive responses

%% repetition suppression

opt = [];
modelClass = [];
rsoa = 10; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

opt.stimContrasts = [.64; .64];
opt.stimDur = 500;
opt.scaling1 = 5e5;
opt.scaling2 = 6e5;
% opt.aAI = 1e2;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 100;
opt.tauS1 = 50;
opt.aAV = 2e2;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

tauE_list = [50 100:100:800];

for ii=1:length(tauE_list)
    % two stim
    opt.stimContrasts = [.64; .64];
    opt.tauE1 = tauE_list(ii);
    opt.tauS1 = 0;
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_E2(:,:,ii) = p.r1;
    
    % one stim
    opt.stimContrasts = [0; .64];
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_E1(:,:,ii) = p.r1;

    % two stim
    opt.stimContrasts = [.64; .64];
    opt.tauE1 = 0;
    opt.tauS1 = tauE_list(ii);
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_S2(:,:,ii) = p.r1;
    
    % one stim
    opt.stimContrasts = [0; .64];
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_S1(:,:,ii) = p.r1;
end


%% backwards masking
