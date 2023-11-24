function analyzeTempWindows

opt = [];
modelClass = [];
rsoa = 4; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2

% H/H, H/L, L/H, L/L
opt.stimContrasts = [.64 .64 .32 .32; ...
                     .64 .32 .64 .32];
% opt.scaling1 = 4e5;
% opt.scaling2 = 5e5;
opt.aAI = 0;

opt.tauE2 = 0;
opt.tauS2 = 0;

opt.display.plotTS = 0; % plot the time series for each simulation

tauE_list = 0:50:1000;
tauS_list = 0:50:1000;

all_taus = combvec(1:length(tauE_list),1:length(tauS_list));
all_perf = nan(length(all_taus),2,4);

parfor ii=1:length(all_taus)
    opt2 = opt;
    
    % loop over taus
    this_tauE = all_taus(1,ii);
    this_tauS = all_taus(2,ii);
    opt2.tauE1 = tauE_list(this_tauE);
    opt2.tauS1 = tauS_list(this_tauS);

    [~,~,perf] = runModel(opt2, modelClass, rsoa, rseq, rcond);
    all_perf(ii,:,:) = perf;
end

all_perf = reshape(all_perf,[length(tauE_list),length(tauS_list),2,4]);

idx=1;
figure
for rr=1:2
  for cc=1:4
    subplot(2,4,idx)
    imagesc(tauE_list,tauS_list,all_perf(:,:,rr,cc));
    %clim([0 1]);
    axis square
    idx = idx+1;
  end
end

%% basic model performance
% high/high
figure
subplot(121)
imagesc(tauE_list,tauS_list,all_perf(:,:,1,1))
clim([0 0.5])
axis square

subplot(122)
imagesc(tauE_list,tauS_list,all_perf(:,:,2,1))
clim([0 0.4])
axis square

% high/low
figure
subplot(121)
imagesc(tauE_list,tauS_list,all_perf(:,:,1,2))
clim([0 0.5])
axis square

subplot(122)
imagesc(tauE_list,tauS_list,all_perf(:,:,2,2))
clim([0 0.16])
axis square

% low/high
figure
subplot(121)
imagesc(tauE_list,tauS_list,all_perf(:,:,1,3))
clim([0 0.22])
axis square

subplot(122)
imagesc(tauE_list,tauS_list,all_perf(:,:,2,3))
clim([0 0.45])
axis square

% low/low
figure
subplot(121)
imagesc(tauE_list,tauS_list,all_perf(:,:,1,4))
clim([0 0.22])
axis square

subplot(122)
imagesc(tauE_list,tauS_list,all_perf(:,:,2,4))
clim([0 0.19])
axis square

%% comparative model performance
% T1 @ high contrast as a function of T2 contrast
figure
imagesc(tauE_list,tauS_list,all_perf(:,:,1,2)-all_perf(:,:,1,1))
clim([0 0.048])
axis square

% T1 @ low contrast as a function of T2 contrast
figure
imagesc(tauE_list,tauS_list,all_perf(:,:,1,4)-all_perf(:,:,1,3))
clim([0 0.025])
axis square

% T2 @ high contrast as a function of T1 contrast
figure
imagesc(tauE_list,tauS_list,all_perf(:,:,2,3)-all_perf(:,:,2,1))
clim([0 0.1])
axis square

% T2 @ low contrast as a function of T1 contrast
figure
imagesc(tauE_list,tauS_list,all_perf(:,:,2,4)-all_perf(:,:,2,2))
clim([0 0.045])
axis square

%% Interactions
% T1
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,1,2)-all_perf(:,:,1,1))-...
    (all_perf(:,:,1,4)-all_perf(:,:,1,3)))
clim([0 0.025])
axis square

% T2
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,2,3)-all_perf(:,:,2,1))-...
    (all_perf(:,:,2,4)-all_perf(:,:,2,2)))
clim([0 0.05])
axis square

% T1 vs T2 - high contrast
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,2,3)-all_perf(:,:,2,1))-...
    (all_perf(:,:,1,2)-all_perf(:,:,1,1)))
clim([-.1 .1])
axis square

% T1 vs T2 - low contrast
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,2,4)-all_perf(:,:,2,2))-...
    (all_perf(:,:,1,4)-all_perf(:,:,1,3)))
clim([-.05 .05])
axis square

%% T1/T2 ratio

% T1 high / T2 high
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,1,1)./all_perf(:,:,2,1)))
colorbar, clim([1 2.1])
axis square

% T1 high / T2 low
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,1,2)./all_perf(:,:,2,2)))
colorbar, clim([1 5])
axis square

% T1 low / T2 high
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,2,3)./all_perf(:,:,1,3)))
colorbar, clim([1 2.2])
axis square

% T1 low / T2 low
figure
imagesc(tauE_list,tauS_list,(all_perf(:,:,1,4)./all_perf(:,:,2,4)))
colorbar, clim([1 1.7])
axis square
