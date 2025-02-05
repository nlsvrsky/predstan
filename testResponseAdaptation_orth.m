% this codes generates model output used for response adaptation analyses
% with orthogonal stimuli in a sequence and was run on the BU SCC

% data output is available on our OSF repository:

ncores = str2num(getenv("NSLOTS"));
pool = parpool(ncores);

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

opt.dt = 2;
opt.T = 7.0*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% generate model output for each parameter combination
tauList = [0 50 100:100:800];
soaList = [0 400:100:1800];

paramList = combvec(tauList,tauList,soaList);

r1_orth = nan(length(paramList),12,opt.nt);

parfor ii=1:length(paramList)
    opt2 = opt;
    opt2.tauE1 = paramList(1,ii);
    opt2.tauS1 = paramList(2,ii);

    thisSOAval = paramList(3,ii);
    if thisSOAval > 0
        thisSOA = thisSOAval;
        opt2.stimContrasts = [.64; .64];
    else
        thisSOA = 1000;
        opt2.stimContrasts = [.64; 0];
    end

    % rseq=3 uses orthogonal stimulus orientations
    [~,p_orth,~] = runModel(opt2, modelClass, thisSOA, 3, rcond);
    r1_orth(ii,:,:) = p_orth.r1;
end

save('output/respAdapt_orth.mat','r1_orth');
