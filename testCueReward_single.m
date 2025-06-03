% adaptation of D-STAN to cue-reward prediction
% (from the response adaptation test)

ncores = 1;%str2num(getenv("NSLOTS"));
%pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 100; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.1;

% this may be needed...
%opt.stimOnset = 1;
%opt.stimDur = 300;

opt.dt = 2;
opt.T = 4.0*1000; % was 7
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output for each parameter combination

% let's just look at one parameter combination
tauE = 100;
tauS = 50;
soa = 2000;
contrast = .90;

%paramList = combvec(tauList,tauList,soaList);
paramList = [tauE; tauS; soa];

r1_iden = nan(12,opt.nt);

opt2 = opt;
opt2.tauE1 = tauE;
opt2.tauS1 = tauS;
opt2.tauE2 = tauE;
opt2.tauS2 = 50*tauS;

thisSOAval = soa;
if thisSOAval > 0
    thisSOA = thisSOAval;
    opt2.stimContrasts = [contrast; contrast];
else
    thisSOA = 1000;
    opt2.stimContrasts = [.64; 0];
end

% rseq=2 uses identical stimulus orientations
[~,p_iden,~] = runModel(opt2, modelClass, thisSOA, 2, rcond);
r1_iden(:,:) = p_iden.r1;
r2_iden(:,:) = p_iden.r2;

save('output/respAdapt_idenSingle.mat','r1_iden');
