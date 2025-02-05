% this codes generates model output used for contrast-dependent suppression
% analyses and was run on the BU SCC

% data output is available on our OSF repository:

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 250; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.1;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output for each parameter combination
tauList = 0:50:1000;
contrList = [.64 .64 .32; .64 .32 .64];

paramList = combvec(tauList,tauList,contrList);

perf_out = nan(2,length(params));
parfor ii=1:length(params)
    opt2 = opt;
    opt2.stimContrasts = paramList(3:4,ii);
    opt2.tauE1 = paramList(1,ii);
    opt2.tauS1 = paramList(2,ii);
    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    perf_out(:,ii) = perf;
end

perf_out = reshape(perf_out,[2 length(tauList) length(tauList) length(contrList)]);

% see Eqn 13
perf_T1 = squeeze((perf_out(1,:,:,2)-perf_out(1,:,:,1))./(perf_out(1,:,:,2)+perf_out(1,:,:,1)));
perf_T2 = squeeze((perf_out(2,:,:,3)-perf_out(2,:,:,1))./(perf_out(2,:,:,3)+perf_out(2,:,:,1)));

save('output/contrDepSupp.mat','perf_out','perf_T1','perf_T2');
