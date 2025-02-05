%% D-STAN results
%  Example code demonstrating each of the main results from our manuscript

addpath('model');

% default model settings
opt = [];
modelClass = [];
rsoa = 2001; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

opt.aAI = 0;
opt.aAV = 0;
opt.sigma1 = 0.1;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

%% subadditivity

opt.tauE1 = 100;
opt.tauS1 = 50;

opt.dt = 2;
opt.T = 4.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

% loop over different stim durations and extract sensory response
stimDur = [30 60 120 240 480];
r1_dur = nan(12,opt.nt,length(stimDur));
for ii=1:length(stimDur)
    opt.stimDur = stimDur(ii);
    [~,p,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);
    r1_dur(:,:,ii) = p.r1;
end

% responses in best neuron
figure
subplot(131)
plot(p.tlist,squeeze(r1_dur(6,:,:)))
xlim([0 2000])

% summed response over time for each duration
subplot(132)
plot(stimDur,squeeze(sum(r1_dur(6,:,:),2)))
xticks(stimDur)
ylim([0 300])

% quantifying subaddtivity
subplot(133)
plot(stimDur(2:end),squeeze(sum(r1_dur(6,:,2:end),2)./sum(r1_dur(6,:,1:end-1))))
ylim([1 2])
xticks(stimDur)

%% response adaptation
%  first we show response adaptation for repeating an identical stimulus

opt.stimDur = 300;

opt.tauE1 = 100;
opt.tauS1 = 50;

soas = [400; 900];
r1_iden = nan(opt.nt,2,2);
for ii=1:length(soas)
    % present both stimuli
    opt.stimContrasts = [.64; .64];

    % rseq=2 presents an identical stimulus
    [~,p,~] = runModel(opt,modelClass,soas(ii),2,rcond);
    r1_iden(:,1,ii) = p.r1(6,:);

    % present only one stimulus
    opt.stimContrasts = [0; .64];
    [~,p,~] = runModel(opt,modelClass,soas(ii),2,rcond);
    r1_iden(:,2,ii) = p.r1(6,:);
end

% response adaptation in T2 at short SOAs
figure
plot(opt.tlist,r1_iden(:,:,1))
xlim([0 2500])
legend({'T1 present','T1 absent'})

% but not long SOAs
figure
plot(opt.tlist,r1_iden(:,:,2))
xlim([0 2500])
legend({'T1 present','T1 absent'})

%% response adaptation in non-identical stimuli
%  similar to before, but now we run a sequence with different orientations

r1_orth = nan(2,opt.nt,2,2); % need responses from 2 neurons
for ii=1:length(soas)
    % present both stimuli
    opt.stimContrasts = [.64; .64];

    % rseq=3 presents orthogonal stimuli
    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_orth(:,:,1,ii) = p.r1([6 12],:);

    % present only one stimulus
    opt.stimContrasts = [0; .64];
    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_orth(:,:,2,ii) = p.r1([6 12],:);
end

% response adaptation in T2 at short SOAs
figure
plot(opt.tlist,r1_orth(:,:,1,1),...
     opt.tlist,r1_orth(:,:,2,1))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

% but not long SOAs
figure
plot(opt.tlist,r1_orth(:,:,1,2),...
     opt.tlist,r1_orth(:,:,2,2))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

%% backward masking

opt.stimDur = 30;

opt.tauE1 = 100;
opt.tauS1 = 50;

soas = [250; 500];
r1_mask = nan(2,opt.nt,2,2);
for ii=1:length(soas)
    % present both stimuli
    opt.stimContrasts = [.64; .64];

    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_mask(:,:,1,ii) = p.r1([6 12],:);

    % present only one stimulus
    opt.stimContrasts = [.64; 0];
    [~,p,~] = runModel(opt,modelClass,soas(ii),3,rcond);
    r1_mask(:,:,2,ii) = p.r1([6 12],:);
end

% backward masking in T1 at short SOAs
figure
plot(opt.tlist,r1_mask(:,:,1,1),...
     opt.tlist,r1_mask(:,:,2,1))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

% but not long SOAs
figure
plot(opt.tlist,r1_mask(:,:,1,2),...
     opt.tlist,r1_mask(:,:,2,2))
xlim([0 2500])
legend({'T1, T1 present','T2, T1 present','T1, T1 absent','T2, T1 absent'})

%% Contrast-dependent suppression

opt.stimDur = 30;
rsoa = 250;

opt.tauE1 = 400;
opt.tauS1 = 100;
opt.scaling1 = 1.25e4;
opt.scaling2 = 1.25e4;

r1_supp = nan(2,opt.nt,4);
d_supp = nan(2,4);

% T1 high, T2 high
opt.stimContrasts = [.64; .64];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,1) = p.r1([6 12],:);
d_supp(:,1) = p.ev;

% T1 high, T2 low
opt.stimContrasts = [.64; .32];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,2) = p.r1([6 12],:);
d_supp(:,2) = p.ev;

% T1 low, T2 high
opt.stimContrasts = [.32; .64];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,3) = p.r1([6 12],:);
d_supp(:,3) = p.ev;

% T1 low, T2 low
opt.stimContrasts = [.32; .32];
[~,p,~] = runModel(opt,modelClass,rsoa,3,rcond);
r1_supp(:,:,4) = p.r1([6 12],:);
d_supp(:,4) = p.ev;


% sensory response as a function of T2 contrast
figure
plot(opt.tlist,r1_supp(1,:,1),'r-',...
     opt.tlist,r1_supp(2,:,1),'b-',...
     opt.tlist,r1_supp(1,:,2),'r--',...
     opt.tlist,r1_supp(2,:,2),'b--')
xlim([0 2500])
legend({'T1, T2 high','T2, T2 high','T1, T2 low','T2, T2 low'})

% sensory response as a function of T1 contrast
figure
plot(opt.tlist,r1_supp(1,:,1),'r-',...
     opt.tlist,r1_supp(2,:,1),'b-',...
     opt.tlist,r1_supp(1,:,3),'r--',...
     opt.tlist,r1_supp(2,:,3),'b--')
xlim([0 2500])
legend({'T1, T1 high','T2, T1 high','T1, T1 low','T2, T1 low'})

% model d' as a function of each contrast level
figure
subplot(121), title('T1')
plot(1:2,d_supp(1,[3 1]),'r-',...
     1:2,d_supp(1,[4 2]),'r--')
axis([[.8 2.2] [0 2.5]])
xticks(1:2)
xticklabels({'T1 low','T1 high'})
legend({'T2 high','T2 low'})

subplot(122), title('T2')
plot(1:2,d_supp(2,[2 1]),'b-',...
     1:2,d_supp(2,[4 3]),'b--')
axis([[.8 2.2] [0 2.5]])
xticks(1:2)
xticklabels({'T2 low','T2 high'})
legend({'T1 high','T1 low'})
