% test stimulus responses

opt = [];
modelClass = [];
rsoa = 2001; % SOA = 250 ms (see runModel)
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
opt.stimDur = 2000;
opt.scaling1 = 5e5;
opt.scaling2 = 6e5;
opt.aAI = 0;
opt.aAV = 0;

opt.eScale = 1;
opt.sScale = 1;

opt.tauE1 = 100;
opt.tauS1 = 50;
% opt.tau1 = 2;
opt.sigma1 = .1;

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
% opt.aAV = 2e2;
% opt.AVWeights = [1 1];
% opt.AVNeutralT1Weight = 0.25;
% opt.distributeVoluntary = 0;

opt.dt = 2;
opt.T = 8.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

tauE_list = [50 100:100:800];
r1_E = nan(12,opt.nt,length(tauE_list));
r1_S = nan(12,opt.nt,length(tauE_list));

parfor ii=1:length(tauE_list)
    opt2 = opt;
    opt2.tauE1 = tauE_list(ii);
    opt2.tauS1 = 0;
    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_E(:,:,ii) = p.r1;
    pE(ii) = p;

    opt2.tauS1 = tauE_list(ii);
    opt2.tauE1 = 0;
    [~,p,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    r1_S(:,:,ii) = p.r1;
    pS(ii) = p;
end

figure
subplot(121)
plot(opt.tlist,squeeze(r1_E(6,:,:)./max(r1_E(6,opt.tlist<2500,:),[],2)))
ylim([0 1.2])

subplot(122)
plot(opt.tlist,squeeze(r1_S(6,:,:)./max(r1_S(6,:,:),[],2)))
ylim([0 1.2])

% time to max
figure
subplot(121)
plot(opt.tlist,squeeze(r1_E(6,:,:)./max(r1_E(6,opt.tlist<2500,:),[],2)))
xlim([490 800])

subplot(122)
plot(opt.tlist,squeeze(r1_S(6,:,:)./max(r1_S(6,:,:),[],2)))
xlim([490 800])

[~,r1_Emax_idx] = max(r1_E(6,opt.tlist<2500,:),[],2);
r1_Emax = squeeze(opt.tlist(r1_Emax_idx));

[~,r1_Smax_idx] = max(r1_S(6,:,:),[],2);
r1_Smax = squeeze(opt.tlist(r1_Smax_idx));

figure
subplot(121)
plot(tauE_list,r1_Emax-500,'.-')
xticks(tauE_list)
xlim([0 850])%, ylim([140 180])

subplot(122)
plot(tauE_list,r1_Smax-500,'.-')
xticks(tauE_list)
xlim([0 850]), ylim([30 60])

% excitatory time to half max
[~,r1_Eh] = min(abs(squeeze(r1_E(6,:,:)./max(r1_E(6,opt.tlist<2500,:),[],2))-0.5));

figure
plot(tauE_list,opt.tlist(r1_Eh)-2500,'.-')
xticks(tauE_list)
xlim([0 850])%, ylim([0 3000])

% suppressive stable activity level
figure
plot(tauE_list,squeeze(r1_S(6,opt.tlist==2500,:)./max(r1_S(6,:,:),[],2)),'.-')
xticks(tauE_list)
xlim([0 850]), ylim([0 0.6])

% time to peak / half max / stable
fig1_ttp = nan(2,9);
fig1_tthm = nan(1,9);
for ii=1:length(tauE_list)
    fig1_ttp(1,ii) = opt.tlist(find(r1_E(6,:,ii)>max(r1_E(6,:,ii))*.99,1))-500;
    fig1_ttp(2,ii) = opt.tlist(find(r1_S(6,:,ii)==max(r1_S(6,:,ii)),1))-500;

    fig1_tthm(ii) = opt.tlist(find(r1_E(6,1250:end,ii)<max(r1_E(6,:,ii))*.5,1));
end

% time to peak
figure
subplot(121)
plot(tauE_list,fig1_ttp(1,:))
axis([0 850 200 1200])

subplot(122)
plot(tauE_list,fig1_ttp(2,:))
axis([0 850 30 60])

% excitatory half max
figure
plot(tauE_list,fig1_tthm)
axis([0 850 0 3000])

% suppressive stable level
figure
plot(tauE_list,squeeze(r1_S(6,1250,:)./max(r1_S(6,:,:),[],2)))
axis([0 850 0 0.6])

% examples of parameter combinations
tauE_list = [100 100 500 500];
tauS_list = [50  500  50 500];

for ii=1:length(tauE_list)
    opt.tauE1 = tauE_list(ii);
    opt.tauS1 = tauS_list(ii);
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_comb(:,:,ii) = p.r1;
end

figure
plot(p.tlist,squeeze(r1_comb(6,:,:)./max(r1_comb(6,:,:),[],2)))
axis([0 8500 0 1.2])


%% subadditive responses
opt = [];
modelClass = [];
rsoa = 1000; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; 0];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;
opt.aAV = 0;

opt.eScale = 1;
opt.sScale = 1;

opt.tauE1 = 100;
opt.tauS1 = 50;
opt.sigma1 = 0.1;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

stimDur = [30 60 120 240 480];
r1_dur = nan(12,opt.nt,length(stimDur));
for ii=1:length(stimDur)
    opt.stimDur = stimDur(ii);
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_dur(:,:,ii) = p.r1;
end

figure
subplot(131)
plot(p.tlist,squeeze(r1_dur(6,:,:)))
xlim([0 2000])

subplot(132)
plot(stimDur,squeeze(sum(r1_dur(6,:,:),2)))
xticks(stimDur)
ylim([0 300])

subplot(133)
plot(stimDur(2:end),squeeze(sum(r1_dur(6,:,2:end),2)./sum(r1_dur(6,:,1:end-1))))
ylim([1 2])
xticks(stimDur)

% effect of time constants
tauE_list = [100 500];
tauS_list = [50  500];

loop_params = combvec(stimDur,tauE_list,tauS_list);

r1_durP = nan(12,opt.nt,length(loop_params));
for ii=1:length(loop_params)
    opt.stimDur = loop_params(1,ii);
    opt.tauE1 = loop_params(2,ii);
    opt.tauS1 = loop_params(3,ii);
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_durP(:,:,ii) = p.r1;
end

r1_durP2 = reshape(r1_durP(6,:,:),length(p.tlist),length(stimDur),length(tauE_list),length(tauS_list));

figure
for ii=1:2, for jj=1:2
    subplot(2,2,(ii-1)*2+jj)
    plot(p.tlist,r1_durP2(:,:,ii,jj))
    title(sprintf("E: %d ms, S: %d ms",tauE_list(ii),tauS_list(jj)))
end, end

r1_durPS = reshape(sum(r1_durP2),5,4);

figure
plot(stimDur(2:end),r1_durPS(2:end,:)./r1_durPS(1:end-1,:))
xticks(stimDur(2:end))
ylim([1 2])

%% repetition suppression

opt = [];
modelClass = [];
rsoa = 1000; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

opt.stimContrasts = [.64; .64];
opt.stimDur = 500;
opt.scaling1 = 5e5;
opt.scaling2 = 6e5;
opt.aAI = 0;
opt.aAV = 0;

opt.eScale = 1;
opt.sScale = 1;

opt.tauE1 = 100;
opt.tauS1 = 50;
opt.sigma1 = 0.1;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

tauE_list = [0 50 100:100:800];

% opt.dt = 2;
% opt.T = 4.1*1000;
% opt.nt = opt.T/opt.dt+1;
% opt.tlist = 0:opt.dt:opt.T;

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
opt = [];
modelClass = [];
rsoa = []; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.stimContrasts = [.64; .64];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.eScale = 1;
opt.sScale = 50;

% opt.tauE1 = 400;
opt.tauS1 = 100;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% effect of excitatory time constant
tauE_list = 50:50:500;
soas = [100:100:1500 2000 2100];
loop_params = combvec(tauE_list,soas);

opt.tauE1 = 100;
rsoa = 100;

opt.stimContrasts = [.64; .64];
[~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);

figure
subplot(121)
plot(p.tlist,p.r1)
xlim([0 2500])
title('T1 and T2')

opt.stimContrasts = [.64; 0];
[~,p,~] = runModel(opt, modelClass, rsoa, rseq, rcond);

subplot(122)
plot(p.tlist,p.r1)
xlim([0 2500])
title('T1 only')
