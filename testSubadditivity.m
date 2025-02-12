% this codes generates model output used for subadditivity analyses
% and was run on the BU SCC

% data output is available on our OSF repository: https://osf.io/qy9pa/

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 1000; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; 0];
opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.1;

opt.dt = 2;
opt.T = 7.0*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output for each parameter combination
tauList = [0 50 100:100:800];
durList = 30 .* 2.^(0:4);

paramList = combvec(tauList,tauList,durList);

r1_sub = nan(length(paramList),12,opt.nt);

parfor ii=1:length(paramList)
    opt2 = opt;
    opt2.tauE1 = paramList(1,ii);
    opt2.tauS1 = paramList(2,ii);
    opt2.stimDur = paramList(3,ii);

    [~,p,~] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_sub(ii,:,:) = p.r1;
end

save('output/respSubadd.mat','r1_sub');
