
opt = [];
modelClass = [];
rsoa = 10; % SOA = 250 ms (see runModel)
rseq = []; % default orientation sequence
rcond = 3; % cueT1, cueT2
soas = [100:50:500 800];

% H/H, H/L, L/H, L/L
% opt.stimContrasts = [.9; .9];
% opt.stimContrasts = [.64 .64 .16 .16; ...
%                      .64 .16 .64 .16];
opt.stimContrasts = [.64 .64 .32 .32; ...
                     .64 .32 .64 .32];
% opt.stimContrasts = [.64; .64];
opt.scaling1 = 3e5;
opt.scaling2 = 4e5;
opt.aAI = 0;

opt.eScale = 1;
opt.sScale = 50;

opt.tauE1 = 900;
opt.tauS1 = 900;
% opt.tau1 = 2;
% opt.sigma1 = .01;

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
% opt.aAV = 1e2;
% opt.AVWeights = [1 1];
% opt.AVNeutralT1Weight = 0.25;
% opt.distributeVoluntary = 1;

opt.display.plotTS = 0; % plot the time series for each simulation
opt.display.plotPerf = 0;

[~,~,perf] = runModel(opt, modelClass, rsoa, rseq, rcond);

if isempty(rcond)
    perf = squeeze(perf);
    
    if size(opt.stimContrasts,2)>1
        sortedPerf = nan(size(perf));
        
        sortedPerf(1,1,:) = perf(1,1,:); % T1 valid
        sortedPerf(1,2,:) = perf(1,3,:); % T1 neutral
        sortedPerf(1,3,:) = perf(1,2,:); % T1 invalid
        
        sortedPerf(2,1,:) = perf(2,2,:); % T2 valid
        sortedPerf(2,2,:) = perf(2,3,:); % T2 neutral
        sortedPerf(2,3,:) = perf(2,1,:); % T2 invalid
        conds = ["T1 high/T2 high","T1 high/T2 low","T1 low/T2 high","T1 low/T2 low"];
        figure
        for ii=1:4
            subplot(1,4,ii)
            plot(1:3, sortedPerf(:,:,ii),'.-','LineWidth',1.5)
            title(conds(ii))
            xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
            ylabel("Performance (d')"), ylim([0 3])
        end
        legend(["T1","T2"])
    
        % compare target performance as a function of non-target
        figure
        subplot(141), hold on, title("T1 high")
        plot(1:3, squeeze(sortedPerf(1,:,[1 2])),'.-','LineWidth',1.5)
        xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
        ylabel("Performance (d')"), ylim([0 3])
        legend("T2 high","T2 low")

        subplot(142), hold on, title("T1 low")
        plot(1:3, squeeze(sortedPerf(1,:,[3 4])),'.-','LineWidth',1.5)
        xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
        ylabel("Performance (d')"), ylim([0 3])
        legend("T2 high","T2 low")
    
        subplot(143), hold on, title("T2 high")
        plot(1:3, squeeze(sortedPerf(2,:,[1 3])),'.-','LineWidth',1.5)
        xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
        ylabel("Performance (d')"), ylim([0 3])
        legend("T1 high","T1 low")
    
        subplot(144), hold on, title("T2 low")
        plot(1:3, squeeze(sortedPerf(2,:,[2 4])),'.-','LineWidth',1.5)
        xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
        ylabel("Performance (d')"), ylim([0 3])
        legend("T1 high","T1 low")
        
        % neutral condition only
        figure
        subplot(1,2,1), hold on, title("T1")
        plot(1:2, squeeze(sortedPerf(1,3,[3 1]))) % T1 with T2 high contrast
        plot(1:2, squeeze(sortedPerf(1,3,[4 2]))) % T1 with T2 low contrast
        xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
        ylim([0 3]), ylabel("Performance (d')")
        legend(["T2 high","T2 low"])
        
        subplot(1,2,2), hold on, title("T2")
        plot(1:2, squeeze(sortedPerf(2,3,[2 1]))) % T2 with T1 high contrast
        plot(1:2, squeeze(sortedPerf(2,3,[4 3]))) % T2 with T1 low contrast
        xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
        ylim([0 3]), ylabel("Performance (d')")
        legend(["T1 high","T1 low"])
    elseif length(opt.stimContrasts)==4
        conds = ["T1 high/T2 high","T1 high/T2 low","T1 low/T2 high","T1 low/T2 low"];
        figure
        for ii=1:4
            subplot(1,4,ii)
            plot(1:3, sortedPerf(:,:,ii),'.-','LineWidth',1.5)
            title(conds(ii))
            xlim([0.5 3.5]), xticks(1:3), xticklabels({'Valid','Neutral','Invalid'})
            ylabel("Performance (d')"), ylim([0 3])
        end
        legend(["T1","T2"])
        
        % neutral condition only
        figure
        subplot(1,2,1), hold on, title("T1")
        plot(1:2, squeeze(sortedPerf(1,3,[3 1]))) % T1 with T2 high contrast
        plot(1:2, squeeze(sortedPerf(1,3,[4 2]))) % T1 with T2 low contrast
        xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
        ylim([0 3]), ylabel("Performance (d')")
        legend(["T2 high","T2 low"])
        
        subplot(1,2,2), hold on, title("T2")
        plot(1:2, squeeze(sortedPerf(2,3,[2 1]))) % T2 with T1 high contrast
        plot(1:2, squeeze(sortedPerf(2,3,[4 3]))) % T2 with T1 low contrast
        xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
        ylim([0 3]), ylabel("Performance (d')")
        legend(["T1 high","T1 low"])
    elseif length(rsoa)>1
        figure
        subplot(121), plot(soas,squeeze(perf(1,:,[1 3 2]))), title("T1")
        xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 3])
        xlabel("SOA (ms)"), ylabel("Performance (d')")
        
        subplot(122), plot(soas,squeeze(perf(2,:,[2 3 1]))), title("T2")
        xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 3])
        xlabel("SOA (ms)"), ylabel("Performance (d')")
        legend(["Valid","Neutral","Invalid"])
    else
        sortedPerf(1,:) = perf(1,[1 3 2]);
        sortedPerf(2,:) = perf(2,[2 3 1]);
        figure
        bar(sortedPerf)
        xticklabels(["T1","T2"])
        ylim([0 3]), ylabel("Performance (d')")
        legend(["Valid","Neutral","Invalid"])
    end
elseif rcond==3 && length(rsoa)==1
    figure
    subplot(1,2,1), hold on, title("T1")
    plot(1:2, perf(1,[3 1]))
    plot(1:2, perf(1,[4 2]))
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    ylim([0 2.5]), ylabel("Performance (d')")
    legend(["T2 high","T2 low"])
    
    subplot(1,2,2), hold on, title("T2")
    plot(1:2, perf(2,[2 1])) % T2 with T1 high contrast
    plot(1:2, perf(2,[4 3])) % T2 with T1 low contrast
    xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
    ylim([0 2.5]), ylabel("Performance (d')")
    legend(["T1 high","T1 low"])
elseif rcond==3 && length(rsoa)>1
    perf = squeeze(perf);
    % low/high performance by SOA
    % for ii=1:length(rsoa)
    %     figure
    %     subplot(1,2,1), hold on, title(sprintf("T1 - %d ms",soas(rsoa(ii))))
    %     plot(1:2, squeeze(perf(1,ii,[3 1])))
    %     plot(1:2, squeeze(perf(1,ii,[4 2])))
    %     xlim([0.8 2.2]), xticks(1:2), xticklabels(["T1 low","T1 high"])
    %     ylim([0 2.5]), ylabel("Performance (d')")
    %     legend(["T2 high","T2 low"])
    % 
    %     subplot(1,2,2), hold on, title(sprintf("T2 - %d ms",soas(rsoa(ii))))
    %     plot(1:2, squeeze(perf(2,ii,[2 1]))) % T2 with T1 high contrast
    %     plot(1:2, squeeze(perf(2,ii,[4 3]))) % T2 with T1 low contrast
    %     xlim([0.8 2.2]), xticks(1:2), xticklabels(["T2 low","T2 high"])
    %     ylim([0 2.5]), ylabel("Performance (d')")
    %     legend(["T1 high","T1 low"])
    % end

    % change in effects across SOAs
    figure
    subplot(121)
    plot(soas(rsoa),squeeze(perf(1,:,[3 4])))
    title("T1 - low")
    xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 2])
    xlabel("SOA (ms)"), ylabel("Performance (d')")
    legend("T2 high","T2 low")

    subplot(122)
    plot(soas(rsoa),squeeze(perf(1,:,[1 2])))
    title("T1 - high")
    xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 2])
    xlabel("SOA (ms)"), ylabel("Performance (d')")

    figure
    subplot(121)
    plot(soas(rsoa),squeeze(perf(2,:,[2 4])))
    title("T2 - low")
    xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 2])
    xlabel("SOA (ms)"), ylabel("Performance (d')")
    legend("T1 high","T2 low")

    subplot(122)
    plot(soas(rsoa),squeeze(perf(2,:,[1 3])))
    title("T2 - high")
    xlim([soas(rsoa(1))-100 soas(rsoa(end))+50]), ylim([0 2])
    xlabel("SOA (ms)"), ylabel("Performance (d')")
end
