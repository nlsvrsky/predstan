%% repetition suppression and backwards masking
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
opt.T = 5.1*1000;
opt.nt = opt.T/opt.dt+1;
opt.tlist = 0:opt.dt:opt.T;

opt.eScale = 1;
opt.sScale = 50;

% opt.tauE1 = 400;
% opt.tauS1 = 100;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

% effect of excitatory time constant
tauE_list = 0:50:500;
soas = [100:100:1500 2000 2100];
loop_params = combvec(tauE_list,soas);

respE = nan(2,opt.nt,length(loop_params));
parfor ii=1:length(loop_params)
    opt2 = opt;
    opt2.tauE1 = loop_params(1,ii);
    opt2.tauS1 = 0;
    this_soa = loop_params(2,ii);
    if this_soa==2000 % repetition suppression baseline
        opt2.stimContrasts = [0; .64];
    elseif this_soa==2100 % backwards masking baseline
        opt2.stimContrasts = [.64; 0];
    else
        opt2.stimContrasts = [.64; .64];
    end

    [~,p,~] = runModel(opt2, modelClass, this_soa, rseq, rcond);
    respE(:,:,ii) = p.r1([6 12],:);
end

respE = reshape(respE,[2,opt.nt,length(tauE_list),length(soas)]);
respE_T1 = squeeze(sum(respE(1,:,:,1:end-2),2) ./ sum(respE(1,:,:,end),2));
respE_T2 = squeeze(sum(respE(2,:,:,1:end-2),2) ./ sum(respE(2,:,:,end-1),2));

figure
plot(soas(1:end-2),respE_T2)
ylim([0 1.05])
title('Repetition suppression')
xlabel('SOA'), ylabel('Normalized response')

figure
plot(soas(1:end-2),respE_T1)
ylim([0 1.05])
title('Backwards masking')
xlabel('SOA'), ylabel('Normalized response')

% effect of suppressive time constant
tauS_list = 0:50:500;
soas = [100:100:1500 2000 2100];
loop_params = combvec(tauS_list,soas);

respS = nan(2,opt.nt,length(loop_params));
parfor ii=1:length(loop_params)
    opt2 = opt;
    opt2.tauE1 = 0;
    opt2.tauS1 = loop_params(1,ii);
    this_soa = loop_params(2,ii);
    if this_soa==2000 % repetition suppression baseline
        opt2.stimContrasts = [0; .64];
    elseif this_soa==2100 % backwards masking baseline
        opt2.stimContrasts = [.64; 0];
    else
        opt2.stimContrasts = [.64; .64];
    end

    [~,p,~] = runModel(opt2, modelClass, this_soa, rseq, rcond);
    respS(:,:,ii) = p.r1([6 12],:);
end

respS = reshape(respS,[2,opt.nt,length(tauE_list),length(soas)]);
respS_T1 = squeeze(sum(respS(1,:,:,1:end-2),2) ./ sum(respS(1,:,:,end),2));
respS_T2 = squeeze(sum(respS(2,:,:,1:end-2),2) ./ sum(respS(2,:,:,end-1),2));

figure
plot(soas(1:end-2),respS_T2)
ylim([0 1.05])
title('Repetition suppression')
xlabel('SOA'), ylabel('Normalized response')

figure
plot(soas(1:end-2),respS_T1)
ylim([0 1.05])
title('Backwards masking')
xlabel('SOA'), ylabel('Normalized response')
