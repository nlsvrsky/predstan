% adaptation of D-STAN to cue-reward prediction

ncores = 1;%str2num(getenv("NSLOTS"));
%pool = parpool(ncores);

addpath('model');

%% model setup
opt = [];
modelClass = [];
rsoa = 100; 
rseq = []; 
rcond = 3; 

opt.aAI = 0;
opt.aAV = 0;
opt.dt = 2;
opt.T = 12000; 
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;
opt.tauE = 100;
opt.tauS = 50;
opt.tauE1 = opt.tauE;
opt.tauS1 = opt.tauS;
opt.tauE2 = opt.tauE;
opt.tauS2 = 150 * opt.tauS;
opt.stimContrasts = [1; 0];

%% task setup
soas = [600 1500 3750 9375] + 250;
predW = 600 ./ (soas/2) + .1; % predictive weight matrix

%% test all four cues separately
cue_results = []; 
for i = 1:length(soas)
    soa = soas(i);
    opt.predW = predW;
    
    [~,p_iden,~] = runModel(opt, modelClass, soa, i, rcond);
    cue_results(end + 1, :) = p_iden.r2;
end

fig = make_linegraph_multiscale(cue_results, p_iden);
saveas(fig, strcat("graph_multiscale.png"));
