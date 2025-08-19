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

% this may be needed...
%opt.stimOnset = 1;
%opt.stimDur = 300;

opt.dt = 2;
opt.T = 12.0*1000; 
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% %% generate model output for each parameter combination

% let's just look at one parameter combination
tauE = 100;
tauS = 50;
soa = 2000;


contrast_results = []; 
for soa = [600 1500 3750 9375] + 250
    scale_cue = 600/(soa/2) + .1;

    paramList = [tauE; tauS; soa];
    
    r1_iden = nan(4,opt.nt);
    r2_iden = nan(1, opt.nt);
    
    opt2 = opt;
    opt2.sigma1 = 1.4;
    opt2.tauE1 = tauE;
    opt2.tauS1 = tauS;
    opt2.tauE2 = tauE;
    opt2.tauS2 = 150*tauS;
    
    thisSOAval = soa;
    if thisSOAval > 0
        thisSOA = thisSOAval;
        opt2.stimContrasts = [1; 1];
    else
        thisSOA = 1000;
        opt2.stimContrasts = [1; 0];
    end
    opt2.scale_cue = scale_cue;
    
    % rseq=2 uses identical stimulus orientations
    [~,p_iden,~] = runModel(opt2, modelClass, thisSOA, 2, rcond);
    r1_iden(:,:) = p_iden.r1; %do we need these right now?
    r2_iden(:,:) = p_iden.r2;
    contrast_results(end + 1, :) = p_iden.r2;
end

fig = make_linegraph_multiscale(contrast_results, p_iden);
saveas(fig, strcat("graph_multiscale.png"));
