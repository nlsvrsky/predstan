function testRandomSeq
% this function generates model responses for reverse correlation analyses
% and was run through the Boston University SCC - adjustments to the code
% are likely necessary to reproduce the analyses

% data output is available on our OSF repository: https://osf.io/qy9pa/

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);
fprintf("Parpool started\n");

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = []; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;
opt.aAV = 0;

opt.sigma1 = 0.1;

opt.dt = 5;
opt.T = 1*3000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.stimMode = 'random'; % generate random stim sequences

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output to random sequences
nSim = 1e4;

stimList = nan(nSim,opt.nt);
d1_out = nan(nSim,opt.nt);
s1_out = nan(nSim,opt.nt);
f1_out = nan(nSim,opt.nt);
r1_out = nan(nSim,opt.nt);

% code was originally run as an array job, with each parameter combination
% executed and output saved separately
paramList = combvec(0:100:900,0:100:900);
paramVal = str2num(getenv("SGE_TASK_ID"));

opt.tauE1 = paramList(1,paramVal);
opt.tauS1 = paramList(2,paramVal);

rng('shuffle');

parfor ii=1:nSim
    [~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);
    stimList(ii,:) = p.stim(1,:);
    d1_out(ii,:) = p.d1(6,:);
    s1_out(ii,:) = p.s1(6,:);
    f1_out(ii,:) = p.f1(6,:);
    r1_out(ii,:) = p.r1(6,:);
end

fprintf('Sims done\n');

out.r1 = r1_out;
out.d1 = d1_out;
out.s1 = s1_out;
out.f1 = f1_out;
out.stimList = stimList;

save(sprintf('output/randomSeq/rand_out_%03d.mat',paramVal),'out');
