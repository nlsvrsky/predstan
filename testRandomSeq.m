
opt = [];
modelClass = [];
rsoa = 10; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

% H/H, H/L, L/H, L/L
opt.stimContrasts = [.64; .64];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;
opt.aAV = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 0;
opt.tauS1 = 0;

opt.dt = 5;
opt.T = 1*1200;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

nSim = 1e4;

stimList = nan(nSim,opt.nt);
d1_out = nan(nSim,opt.nt);
s1_out = nan(nSim,opt.nt);
f1_out = nan(nSim,opt.nt);
r1_out = nan(nSim,opt.nt);

paramList = combvec(1:nSim,[400; 100]);%[0 100 400; 0 400 100]);

parfor ii=1:length(paramList)
    opt2 = opt;
    opt2.tauE1 = paramList(2,ii);
    opt2.tauS1 = paramList(3,ii);
    [~,p,~] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    stimList(ii,:) = p.stim(1,:);
    d1_out(ii,:) = p.d1(6,:);
    s1_out(ii,:) = p.s1(6,:);
    f1_out(ii,:) = p.f1(6,:);
    r1_out(ii,:) = p.r1(6,:);
end
